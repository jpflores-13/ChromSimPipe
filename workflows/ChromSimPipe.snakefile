#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# ChromSimPipe — cohesin loop-extrusion simulation workflow
#
# Run with:  sbatch ChromSimPipe.sh
#
# All inputs come from the config file (config/ChromSimConfig.yaml) and the
# samplesheet it names. configs/parameters.py parses both; this file only
# wires the rules together.
#
#   convertHic (per Hi-C sample) ─────────────────────────────┐
#   getMotif ─► orientCTCF (per CTCF sample) ─► validateCTCF  │
#                    ─► simulateShard (condition × rep × shard, GPU)
#                    ─► mergeShards (condition × rep)         │
#                    ─► analyzeCondition (condition) ◄────────┘
#                    ─► benchmarkSummary
# ==============================================================================

import os
import sys

##### Load config and sample sheets #####
## `--configfile other.yaml` replaces the default config entirely (it must be
## a complete config); configs/parameters.py reads the same file, here and in
## every job, via $CHROMSIM_CONFIG.
DEFAULT_CONFIG = "config/ChromSimConfig.yaml"
if not workflow.overwrite_configfiles:
    configfile: DEFAULT_CONFIG
CONFIGFILE = os.path.abspath((list(workflow.overwrite_configfiles) or [DEFAULT_CONFIG])[-1])
os.environ["CHROMSIM_CONFIG"] = CONFIGFILE

sys.path.insert(0, os.getcwd())
from configs.parameters import (
    SIMULATION_CONDITIONS, CTCF_TYPES, CTCF_PEAKS, CTCF_BEDS,
    HIC_FILES, MCOOL_FILES, PATHS, SAMPLES, RUN_NAME,
    REGION, N_REPLICATES, N_SHARDS, ANALYSIS, LOCUS_TAG, GENOME_ASSEMBLY,
)


def rel(path):
    """Paths inside the working directory are written relative to it."""
    if not os.path.isabs(path):
        return path
    r = os.path.relpath(path)
    return path if r.startswith("..") else r


def conda_env(key):
    """Resolve a condaEnvs entry: an env file, an existing env name, or an env path."""
    env = os.path.expanduser(str(config["condaEnvs"][key]))
    if env.endswith((".yml", ".yaml")) and not os.path.isabs(env):
        return os.path.abspath(env)   # Snakemake resolves env files from workflows/ otherwise
    return env


COND_LOOKUP = {c["name"]: c for c in SIMULATION_CONDITIONS}
COND_NAMES  = list(COND_LOOKUP)

OUT        = rel(os.path.dirname(PATHS["mcool"]))
MCOOL_DIR  = rel(PATHS["mcool"])
CTCF_DIR   = rel(PATHS["ctcf"])
MOTIF_DIR  = rel(PATHS["motifs"])
SIM_DIR    = rel(PATHS["sims"])
ANA_DIR    = rel(PATHS["analysis"])
LOG_DIR    = rel(PATHS["logs"])
BENCH_DIR  = rel(PATHS["benchmarks"])
DONE_DIR   = rel(PATHS["done"])

ENV_SIM  = conda_env("sim")
ENV_CTCF = conda_env("ctcf")

MOTIF = config["ctcf_motif"]
MOTIF_MEME = (os.path.expanduser(MOTIF["meme"]) if MOTIF.get("meme")
              else f"{MOTIF_DIR}/{MOTIF['jaspar_id']}.meme")

wildcard_constraints:
    rep   = r"\d+",
    shard = r"\d+"

## Write the parsed samplesheet (with sample names) next to the outputs.
## Only rewritten when it changes: every cluster job re-parses this file.
_sheet = f"{OUT}/{RUN_NAME}_ChromSimSamplesheet.txt"
_cols  = list(next(iter(SAMPLES.values()))["metadata"])
_text  = "\t".join(["sampleName"] + _cols) + "\n" + "".join(
    "\t".join([name] + [s["metadata"][c] for c in _cols]) + "\n"
    for name, s in SAMPLES.items()
)
if not os.path.exists(_sheet) or open(_sheet).read() != _text:
    os.makedirs(OUT, exist_ok=True)
    with open(_sheet, "w") as _fh:
        _fh.write(_text)


def compare_mcool(wildcards):
    hic = COND_LOOKUP[wildcards.condition]["compare_hic"]
    return [rel(MCOOL_FILES[hic])] if hic else []


## Define actions on success
onsuccess:
    print("ChromSimPipe completed successfully! Wahoo!")


##### Define rules #####
rule all:
    input:
        expand(f"{DONE_DIR}/analyzeCondition_{{condition}}.done", condition=COND_NAMES),
        f"{OUT}/{RUN_NAME}_benchmarking.tsv"


rule convertHic:
    input:
        hic = lambda wildcards: HIC_FILES[wildcards.sampleName]
    output:
        cool = f"{MCOOL_DIR}/{{sampleName}}.cool"
    log:
        err = f"{LOG_DIR}/convertHic_{{sampleName}}.err",
        out = f"{LOG_DIR}/convertHic_{{sampleName}}.out"
    benchmark:
        f"{BENCH_DIR}/convertHic_{{sampleName}}.tsv"
    params:
        resolution = config["resolution"]
    conda: ENV_SIM
    shell:
        """
        hic2cool convert {input.hic} {output.cool} -r {params.resolution} 1> {log.out} 2> {log.err}
        """


if not MOTIF.get("meme"):
    rule getMotif:
        output:
            meme = MOTIF_MEME
        log:
            err = f"{LOG_DIR}/getMotif_{MOTIF['jaspar_id']}.err",
            out = f"{LOG_DIR}/getMotif_{MOTIF['jaspar_id']}.out"
        benchmark:
            f"{BENCH_DIR}/getMotif_{MOTIF['jaspar_id']}.tsv"
        params:
            url = MOTIF["jaspar_url"].format(jaspar_id=MOTIF["jaspar_id"])
        conda: ENV_CTCF
        shell:
            """
            wget -O {output.meme} "{params.url}" 1> {log.out} 2> {log.err}
            grep -q "^MOTIF" {output.meme} || {{ echo "Not a MEME file: {params.url}" >> {log.err}; exit 1; }}
            """


rule orientCTCF:
    input:
        peaks  = lambda wildcards: CTCF_PEAKS[wildcards.ctcfType],
        genome = config["genome"],
        motif  = MOTIF_MEME
    output:
        bed = f"{CTCF_DIR}/ctcf_oriented_{GENOME_ASSEMBLY}_{{ctcfType}}_{LOCUS_TAG}.bed"
    log:
        err = f"{LOG_DIR}/orientCTCF_{{ctcfType}}.err",
        out = f"{LOG_DIR}/orientCTCF_{{ctcfType}}.out"
    benchmark:
        f"{BENCH_DIR}/orientCTCF_{{ctcfType}}.tsv"
    params:
        region     = REGION,
        assembly   = GENOME_ASSEMBLY,
        thresh     = MOTIF["fimo_thresh"],
        flank      = MOTIF["flank_bp"],
        unoriented = MOTIF["unoriented"],
        resolution = config["resolution"]
    conda: ENV_CTCF
    shell:
        """
        python scripts/extract_ctcf_sites.py \
            --bed {input.peaks} \
            --region {params.region} \
            --genome {input.genome} \
            --motif {input.motif} \
            --assembly {params.assembly} \
            --fimo-thresh {params.thresh} \
            --flank {params.flank} \
            --unoriented {params.unoriented} \
            --resolution {params.resolution} \
            --output {output.bed} 1> {log.out} 2> {log.err}
        """


rule validateCTCF:
    input:
        beds = [rel(CTCF_BEDS[t]) for t in CTCF_TYPES]
    output:
        done = f"{DONE_DIR}/validateCTCF.done"
    log:
        err = f"{LOG_DIR}/validateCTCF.err",
        out = f"{LOG_DIR}/validateCTCF.out"
    benchmark:
        f"{BENCH_DIR}/validateCTCF.tsv"
    params:
        configfile = CONFIGFILE
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} python scripts/validate_ctcf.py 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule simulateShard:
    input:
        validated = rules.validateCTCF.output.done,
        bed       = lambda wildcards: rel(CTCF_BEDS[COND_LOOKUP[wildcards.condition]["ctcf_type"]])
    output:
        done = f"{DONE_DIR}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.done"
    log:
        err = f"{LOG_DIR}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.err",
        out = f"{LOG_DIR}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.out"
    benchmark:
        f"{BENCH_DIR}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = SIM_DIR,
        n_shards   = N_SHARDS
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} python scripts/run_simulation_shard.py \
            --condition {wildcards.condition} \
            --replicate {wildcards.rep} \
            --shard-index {wildcards.shard} \
            --n-shards {params.n_shards} \
            --gpu 0 \
            --output {params.sim_dir} 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule mergeShards:
    input:
        lambda wildcards: expand(rules.simulateShard.output.done,
                                 condition=wildcards.condition, rep=wildcards.rep,
                                 shard=range(N_SHARDS))
    output:
        done = f"{DONE_DIR}/mergeShards_{{condition}}_rep{{rep}}.done"
    log:
        err = f"{LOG_DIR}/mergeShards_{{condition}}_rep{{rep}}.err",
        out = f"{LOG_DIR}/mergeShards_{{condition}}_rep{{rep}}.out"
    benchmark:
        f"{BENCH_DIR}/mergeShards_{{condition}}_rep{{rep}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = SIM_DIR,
        n_shards   = N_SHARDS
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} python scripts/merge_shards.py \
            --condition {wildcards.condition} \
            --replicate {wildcards.rep} \
            --n-shards {params.n_shards} \
            --results-dir {params.sim_dir} 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule analyzeCondition:
    input:
        merged = lambda wildcards: expand(rules.mergeShards.output.done,
                                          condition=wildcards.condition,
                                          rep=range(N_REPLICATES)),
        mcool  = compare_mcool
    output:
        done = f"{DONE_DIR}/analyzeCondition_{{condition}}.done"
    log:
        err = f"{LOG_DIR}/analyzeCondition_{{condition}}.err",
        out = f"{LOG_DIR}/analyzeCondition_{{condition}}.out"
    benchmark:
        f"{BENCH_DIR}/analyzeCondition_{{condition}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = SIM_DIR,
        out_dir    = lambda wildcards: f"{ANA_DIR}/{wildcards.condition}",
        mcool_dir  = MCOOL_DIR,
        n_jobs     = ANALYSIS["n_jobs"],
        msd_flag   = "" if ANALYSIS["msd"] else "--no-msd"
    conda: ENV_SIM
    shell:
        """
        mkdir -p {params.out_dir}
        CHROMSIM_CONFIG={params.configfile} python scripts/run_analysis_all.py \
            --results-dir {params.sim_dir} \
            --output-dir {params.out_dir} \
            --condition {wildcards.condition} \
            --hic-dir {params.mcool_dir} \
            --n-jobs {params.n_jobs} \
            --skip-existing \
            {params.msd_flag} 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule benchmarkSummary:
    input:
        expand(rules.analyzeCondition.output.done, condition=COND_NAMES)
    output:
        f"{OUT}/{RUN_NAME}_benchmarking.tsv"
    log:
        err = f"{LOG_DIR}/benchmarkSummary.err"
    params:
        bench_dir = BENCH_DIR
    conda: ENV_SIM
    shell:
        """
        python workflows/utils/benchmarking.py {params.bench_dir} {output} 2> {log.err}
        """
