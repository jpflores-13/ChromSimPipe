# ==============================================================================
# ChromSimPipe Snakemake workflow
# Cohesin loop-extrusion simulation pipeline
#
# Usage:  sbatch ChromSimPipe.sh
#         (or see ChromSimPipe.sh for the full snakemake command)
#
# Inputs (set in config/ChromSimConfig.yaml):
#   hic_control / hic_sorbitol   — .hic contact maps (one per condition)
#   genome                        — genome FASTA for FIMO orientation
#   samplesheet                   — ChromSimSamplesheet.txt (CTCF peaks per condition)
#
# Locus + conditions are read from configs/parameters.py (single source of truth).
# Change the locus section in ChromSimConfig.yaml to switch genomic regions.
# ==============================================================================

import sys, os
import pandas as pd

sys.path.insert(0, os.getcwd())
from configs.parameters import (
    SIMULATION_CONDITIONS, LOCI, SIM_RUN, ACTIVE_LOCUS
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
        hic = lambda wc: config["hic_" + wc.hic_type]
    output:
        mcool = "data/mcool/{hic_type}.cool"
    resources:
        runtime        = 120,
        cpus_per_task  = 4,
        mem_mb         = 32000,
    shell:
        """
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p data/mcool
        conda run -n cohesin_sim hic2cool convert \\
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
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p data/ctcf_beds
        genome_arg=""
        [ -f "{params.genome}" ] && genome_arg="--genome {params.genome}"
        conda run -n ctcf_extraction python scripts/extract_ctcf_sites_hg38.py \\
            --source bed \\
            --bed {input.peaks} \\
            --region {params.region} \\
            $genome_arg \\
            --output {output.bed}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — Validate CTCF BED files (runs once, after both conditions extracted)
# ══════════════════════════════════════════════════════════════════════════════
rule validate_ctcf:
    input:
        control  = _ctcf_bed("control"),
        sorbitol = _ctcf_bed("sorbitol"),
    output:
        done = SENTINEL + "/ctcf_validated.done"
    resources:
        runtime        = 15,
        cpus_per_task  = 1,
        mem_mb         = 4000,
    shell:
        """
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p $(dirname {output.done})
        conda run -n cohesin_sim python scripts/validate_ctcf.py
        touch {output.done}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — GPU simulation: one SLURM job per (condition × replicate × shard)
#   60 GPU jobs total: 5 conditions × 3 replicates × 4 shards
# ══════════════════════════════════════════════════════════════════════════════
rule simulate_shard:
    input:
        ctcf_done      = SENTINEL + "/ctcf_validated.done",
        control_mcool  = "data/mcool/control.cool",
        sorbitol_mcool = "data/mcool/sorbitol.cool",
    output:
        done = SENTINEL + "/sims/{condition}_rep{rep}_shard{shard}.done"
    params:
        results_dir = RESULTS_DIR,
        n_shards    = N_SHARDS,
    resources:
        runtime         = 720,
        slurm_partition = "gpu",
        cpus_per_task   = 4,
        mem_mb          = 64000,
        slurm_extra     = "'--gpus=1'",
    shell:
        """
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p {params.results_dir}
        mkdir -p $(dirname {output.done})
        conda run -n cohesin_sim python scripts/run_simulation_shard.py \\
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
#   15 merge jobs: 5 conditions × 3 replicates
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
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p $(dirname {output.done})
        conda run -n cohesin_sim python scripts/merge_shards.py \\
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
        control_mcool  = "data/mcool/control.cool",
        sorbitol_mcool = "data/mcool/sorbitol.cool",
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
        module load anaconda/2024.02 2>/dev/null || true
        mkdir -p {params.output_dir}
        mkdir -p $(dirname {output.done})
        conda run -n cohesin_sim python scripts/run_analysis_all.py \\
            --results-dir {params.results_dir} \\
            --output-dir {params.output_dir} \\
            --condition {wildcards.condition} \\
            --hic-dir {params.hic_dir} \\
            --mcool-mesc {input.control_mcool} \\
            --mcool-neuron {input.sorbitol_mcool} \\
            --n-jobs {params.n_jobs} \\
            --skip-existing \\
            --no-msd
        touch {output.done}
        """
