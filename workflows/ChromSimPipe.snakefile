# ==============================================================================
# ChromSimPipe Snakemake workflow
# Cohesin loop-extrusion simulation pipeline
#
# Usage:  sbatch ChromSimPipe.sh
#         (or see ChromSimPipe.sh for the full snakemake command)
#
# Inputs (set in config/ChromSimConfig.yaml):
#   hic                           — {hic_type: .hic path}, one per Hi-C map
#   genome                        — genome FASTA for FIMO orientation
#   samplesheet                   — ChromSimSamplesheet.txt (CTCF peaks per ctcf_type)
#
# Locus + conditions are read from configs/parameters.py (single source of truth).
# Change the locus section in ChromSimConfig.yaml to switch genomic regions.
# ==============================================================================

import sys, os
import pandas as pd

sys.path.insert(0, os.getcwd())
from configs.parameters import (
    SIMULATION_CONDITIONS, LOCI, SIM_RUN, ACTIVE_LOCUS,
    CTCF_TYPES, CTCF_BEDS, HIC_FILES, MCOOL_FILES,
)

configfile: "config/ChromSimConfig.yaml"

# ── Pipeline-level constants ───────────────────────────────────────────────────
COND_LOOKUP  = {c["name"]: c for c in SIMULATION_CONDITIONS}
COND_NAMES   = [c["name"] for c in SIMULATION_CONDITIONS]

LOCUS      = ACTIVE_LOCUS
LOCUS_INFO = LOCI[LOCUS]
CHROM      = LOCUS_INFO["chrom"]
RSTART     = LOCUS_INFO["region_start"]
REND       = LOCUS_INFO["region_end"]
REGION     = f"{CHROM}:{RSTART}-{REND}"

N_REP        = config.get("n_replicates", SIM_RUN["n_replicates"])
N_SHARDS     = config.get("n_shards", 4)
RESULTS_DIR  = "results/polychrom_3d"
SENTINEL     = "results/.sentinels"

# Read samplesheet (CTCF_Type → CTCF_Peaks_Path)
_samples = pd.read_table(config["samplesheet"]).set_index("CTCF_Type", drop=False)

def peaks_path(ctcf_type):
    return _samples.loc[ctcf_type, "CTCF_Peaks_Path"]

_ASSEMBLY  = config.get("genome_assembly", "hg38")
_locus_tag = f"{CHROM}_{RSTART}_{REND}"

def _ctcf_bed(ctcf_type):
    return f"data/ctcf_beds/ctcf_oriented_{_ASSEMBLY}_{ctcf_type}_{_locus_tag}.bed"

assert all(CTCF_BEDS[t] == _ctcf_bed(t) for t in CTCF_TYPES), \
    "CTCF BED naming differs between the Snakefile and configs/parameters.py"

_missing_peaks = [t for t in CTCF_TYPES if t not in _samples.index]
if _missing_peaks:
    raise ValueError(f"No samplesheet row for ctcf_type(s) {_missing_peaks} "
                     f"in {config['samplesheet']}")

def _compare_mcool(wc):
    hic = COND_LOOKUP[wc.condition]["compare_hic"]
    return [MCOOL_FILES[hic]] if hic else []


# ══════════════════════════════════════════════════════════════════════════════
# rule all — final target
# ══════════════════════════════════════════════════════════════════════════════
rule all:
    input:
        expand(SENTINEL + "/analyzed/{condition}.done", condition=COND_NAMES)


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 — Convert .hic → .mcool (one job per Hi-C map)
# ══════════════════════════════════════════════════════════════════════════════
rule convert_hic:
    input:
        hic = lambda wc: HIC_FILES[wc.hic_type]
    output:
        mcool = "data/mcool/{hic_type}.cool"
    resources:
        runtime        = 120,
        cpus_per_task  = 4,
        mem_mb         = 32000,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p data/mcool
        conda run -p ~/.local/share/mamba/envs/cohesin_sim hic2cool convert \\
            {input.hic} {output.mcool} -r {config[resolution]}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — Extract oriented CTCF sites via FIMO (one job per CTCF condition)
# ══════════════════════════════════════════════════════════════════════════════
rule extract_ctcf:
    input:
        peaks = lambda wc: peaks_path(wc.ctcf_type)
    output:
        bed = f"data/ctcf_beds/ctcf_oriented_{_ASSEMBLY}_{{ctcf_type}}_{_locus_tag}.bed"
    params:
        region = REGION,
        genome = config.get("genome", ""),
    resources:
        runtime        = 60,
        cpus_per_task  = 2,
        mem_mb         = 16000,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p data/ctcf_beds
        genome_arg=""
        [ -f "{params.genome}" ] && genome_arg="--genome {params.genome}"
        conda run -p ~/.local/share/mamba/envs/ctcf_extraction python scripts/extract_ctcf_sites_hg38.py \\
            --source bed \\
            --bed {input.peaks} \\
            --region {params.region} \\
            $genome_arg \\
            --output {output.bed}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — Validate CTCF BED files (runs once, after every ctcf_type is extracted)
# ══════════════════════════════════════════════════════════════════════════════
rule validate_ctcf:
    input:
        beds = [_ctcf_bed(t) for t in CTCF_TYPES],
    output:
        done = SENTINEL + "/ctcf_validated.done"
    resources:
        runtime        = 15,
        cpus_per_task  = 1,
        mem_mb         = 4000,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p $(dirname {output.done})
        conda run -p ~/.local/share/mamba/envs/cohesin_sim python scripts/validate_ctcf.py
        touch {output.done}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — GPU simulation: one SLURM job per (condition × replicate × shard)
#   (conditions × n_replicates × n_shards) GPU jobs
# ══════════════════════════════════════════════════════════════════════════════
rule simulate_shard:
    input:
        ctcf_done = SENTINEL + "/ctcf_validated.done",
        ctcf_bed  = lambda wc: _ctcf_bed(COND_LOOKUP[wc.condition]["ctcf_type"]),
    output:
        done = SENTINEL + "/sims/{condition}_rep{rep}_shard{shard}.done"
    params:
        results_dir = RESULTS_DIR,
        n_shards    = N_SHARDS,
    resources:
        runtime         = 720,
        slurm_partition = "gpu",
        qos             = "gpu_access",
        cpus_per_task   = 4,
        mem_mb          = 60000,
        mem_mb_per_cpu  = 0,
        gpu             = 1,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p {params.results_dir}
        mkdir -p $(dirname {output.done})
        conda run -p ~/.local/share/mamba/envs/cohesin_sim python scripts/run_simulation_shard.py \\
            --condition {wildcards.condition} \\
            --replicate {wildcards.rep} \\
            --shard-index {wildcards.shard} \\
            --n-shards {params.n_shards} \\
            --gpu 0 \\
            --output {params.results_dir}
        touch {output.done}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 — Merge shards: one CPU job per (condition × replicate)
#   (conditions × n_replicates) merge jobs
# ══════════════════════════════════════════════════════════════════════════════
rule merge_shards:
    input:
        lambda wc: expand(
            SENTINEL + "/sims/{condition}_rep{rep}_shard{shard}.done",
            condition=wc.condition,
            rep=wc.rep,
            shard=range(N_SHARDS)
        )
    output:
        done = SENTINEL + "/merged/{condition}_rep{rep}.done"
    params:
        results_dir = RESULTS_DIR,
        n_shards    = N_SHARDS,
    resources:
        runtime        = 120,
        cpus_per_task  = 4,
        mem_mb         = 64000,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p $(dirname {output.done})
        conda run -p ~/.local/share/mamba/envs/cohesin_sim python scripts/merge_shards.py \\
            --condition {wildcards.condition} \\
            --replicate {wildcards.rep} \\
            --n-shards {params.n_shards} \\
            --results-dir {params.results_dir}
        touch {output.done}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 6 — Analysis: contact maps, P(s) curves, MSD, figures (one job per cond)
# ══════════════════════════════════════════════════════════════════════════════
rule analyze_all:
    input:
        merged = lambda wc: expand(
            SENTINEL + "/merged/{condition}_rep{rep}.done",
            condition=wc.condition,
            rep=range(N_REP)
        ),
        mcool = _compare_mcool,
    output:
        done = SENTINEL + "/analyzed/{condition}.done"
    params:
        results_dir = RESULTS_DIR,
        output_dir  = lambda wc: f"results/analysis/{wc.condition}",
        hic_dir     = "data/mcool",
        n_jobs      = 4,
    resources:
        runtime        = 360,
        cpus_per_task  = 4,
        mem_mb         = 64000,
    shell:
        """
        source /nas/longleaf/rhel9/apps/anaconda/2024.02/etc/profile.d/conda.sh
        mkdir -p {params.output_dir}
        mkdir -p $(dirname {output.done})
        conda run -p ~/.local/share/mamba/envs/cohesin_sim python scripts/run_analysis_all.py \\
            --results-dir {params.results_dir} \\
            --output-dir {params.output_dir} \\
            --condition {wildcards.condition} \\
            --hic-dir {params.hic_dir} \\
            --n-jobs {params.n_jobs} \\
            --skip-existing \\
            --no-msd
        touch {output.done}
        """
