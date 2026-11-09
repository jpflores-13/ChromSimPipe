#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# ChromSimPipe - cohesin loop-extrusion simulation workflow (polysim model)
#
# Run with:  sbatch ChromSimPipe.sh
#
# All inputs come from the config file (config/ChromSimConfig.yaml) and the
# samplesheet it names. configs/parameters.py parses both; this file only
# wires the rules together. Every locus in the config's `loci:` list is
# simulated and analysed separately, in <outputDir>/loci/<locus>/.
#
#   Rules, in order:
#     convertHic        per Hi-C sample (shared by all loci)
#     getMotif          once
#     preparePolysim    once (builds polysim's Cython LEF model)
#     orientCTCF        per locus x CTCF sample
#     validateCTCF      per locus
#     simulateShard     per locus x condition x replicate x shard (GPU)
#     mergeShards       per locus x condition x replicate
#     analyzeCondition  per locus x condition (needs convertHic output)
#     benchmarkSummary  once
# ==============================================================================

import os
import sys

##### Load config and sample sheets #####
## `--configfile other.yaml` replaces the default config entirely (it must be
## a complete config); configs/parameters.py reads the same file, here and in
## every job, via $CHROMSIM_CONFIG. Each per-locus job also gets $CHROMSIM_LOCUS.
DEFAULT_CONFIG = "config/ChromSimConfig.yaml"
if not workflow.overwrite_configfiles:
    configfile: DEFAULT_CONFIG
CONFIGFILE = os.path.abspath((list(workflow.overwrite_configfiles) or [DEFAULT_CONFIG])[-1])
os.environ["CHROMSIM_CONFIG"] = CONFIGFILE

sys.path.insert(0, os.getcwd())
from configs.parameters import (
    SIMULATION_CONDITIONS, CTCF_TYPES, CTCF_PEAKS, HIC_FILES, MCOOL_FILES,
    SAMPLES, RUN_NAME, N_REPLICATES, N_SHARDS, ANALYSIS, CALIBRATION, GENOME_ASSEMBLY, GENOME,
    OUTPUT_DIR, LOCI, LOCUS_NAMES, locus_paths, ctcf_bed_path, condition_peak_files,
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

OUT        = rel(OUTPUT_DIR)
MCOOL_DIR  = rel(locus_paths(LOCUS_NAMES[0])["mcool"])
MOTIF_DIR  = rel(locus_paths(LOCUS_NAMES[0])["motifs"])
LOG_DIR    = f"{OUT}/logs"
BENCH_DIR  = f"{OUT}/benchmarks"


def L(key):
    """A per-locus folder as a pattern with a {locus} wildcard."""
    return f"{OUT}/loci/{{locus}}/{key}"


def LP(wildcards, key):
    """A per-locus folder for a job's locus."""
    return rel(locus_paths(wildcards.locus)[key])


ENV_SIM  = conda_env("sim")
ENV_CTCF = conda_env("ctcf")

MOTIF = config["ctcf_motif"]
MOTIF_MEME = (os.path.expanduser(MOTIF["meme"]) if MOTIF.get("meme")
              else f"{MOTIF_DIR}/{MOTIF['jaspar_id']}.meme")

## Rules call "$CONDA_PREFIX/bin/python" (the rule's conda env) rather than a
## bare `python`: in cluster jobs another interpreter (e.g. the launcher's
## venv) can come first on PATH and lacks the env's packages.

wildcard_constraints:
    locus = r"[^/]+",
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


## Email notifications (notify: in the config; $CHROMSIM_EMAIL overrides the address)
NOTIFY = config.get("notify") or {}
EMAIL  = os.environ.get("CHROMSIM_EMAIL") or NOTIFY.get("email") or ""


def notify(event, body):
    """Send one message for a workflow event, if an address and the event are configured."""
    if not EMAIL or event not in NOTIFY.get("events", []):
        return
    import shlex, socket, subprocess
    subject = f"ChromSimPipe {RUN_NAME}: {event}"
    body = (f"Run:     {RUN_NAME}\nConfig:  {CONFIGFILE}\nOutput:  {os.path.abspath(OUT)}\n"
            f"Host:    {socket.gethostname()}\n\n{body}\n")
    cmd = NOTIFY["command"].format(subject=shlex.quote(subject), email=shlex.quote(EMAIL))
    try:
        subprocess.run(cmd, shell=True, input=body, text=True, timeout=60, check=True)
    except Exception as e:  # a failed email must never fail the workflow
        print(f"Could not send '{event}' email to {EMAIL}: {e}", file=sys.stderr)


def log_tail(path, n=40):
    try:
        with open(path) as fh:
            return "".join(fh.readlines()[-n:])
    except (OSError, TypeError):
        return "(Snakemake log unavailable)"


## Define actions on start, success, and error
onstart:
    notify("start", f"Started. Loci: {', '.join(LOCUS_NAMES)}. Conditions: {', '.join(COND_NAMES)}")

onsuccess:
    print("ChromSimPipe completed successfully! Wahoo!")
    notify("success", f"Finished successfully.\nResults: {os.path.abspath(OUT)}/loci/<locus>/analysis\n"
                      f"Timings: {os.path.abspath(OUT)}/{RUN_NAME}_benchmarking.tsv")

onerror:
    notify("error", f"A job failed. Per-job logs: {os.path.abspath(LOG_DIR)}\n"
                    f"Snakemake log: {log}\n\nLast lines:\n{log_tail(log)}")


##### Define rules #####
rule all:
    input:
        expand(f"{L('done')}/analyzeCondition_{{condition}}.done",
               locus=LOCUS_NAMES, condition=COND_NAMES),
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
        "$CONDA_PREFIX/bin/python" scripts/prepare_hic.py --input {input.hic} --output {output.cool} \
            --resolution {params.resolution} 1> {log.out} 2> {log.err}
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


## polysim builds its Cython LEF model (LEF_Dynamics.pyx) with pyximport on
## first import. Do that once here, so parallel shards load the built module
## instead of compiling it at the same time.
rule preparePolysim:
    output:
        done = f"{OUT}/done/preparePolysim.done"
    log:
        err = f"{LOG_DIR}/preparePolysim.err"
    conda: ENV_SIM
    shell:
        """
        "$CONDA_PREFIX/bin/python" -c "import polysim.extrusion, polysim.sim3d" 2> {log.err}
        touch {output.done}
        """


rule orientCTCF:
    input:
        peaks  = lambda wildcards: CTCF_PEAKS[wildcards.ctcfType],
        genome = rel(GENOME),
        motif  = MOTIF_MEME
    output:
        bed = f"{L('ctcf')}/ctcf_oriented_{GENOME_ASSEMBLY}_{{ctcfType}}.bed"
    log:
        err = f"{LOG_DIR}/{{locus}}/orientCTCF_{{ctcfType}}.err",
        out = f"{LOG_DIR}/{{locus}}/orientCTCF_{{ctcfType}}.out"
    benchmark:
        f"{BENCH_DIR}/orientCTCF_{{locus}}_{{ctcfType}}.tsv"
    params:
        region     = lambda wildcards: LOCI[wildcards.locus]["region"],
        assembly   = GENOME_ASSEMBLY,
        thresh     = MOTIF["fimo_thresh"],
        flank      = MOTIF["flank_bp"],
        unoriented = MOTIF["unoriented"],
        resolution = config["resolution"]
    conda: ENV_CTCF
    shell:
        """
        "$CONDA_PREFIX/bin/python" scripts/extract_ctcf_sites.py \
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
        beds = lambda wildcards: [rel(ctcf_bed_path(wildcards.locus, t)) for t in CTCF_TYPES]
    output:
        done = f"{L('done')}/validateCTCF.done"
    log:
        err = f"{LOG_DIR}/{{locus}}/validateCTCF.err",
        out = f"{LOG_DIR}/{{locus}}/validateCTCF.out"
    benchmark:
        f"{BENCH_DIR}/validateCTCF_{{locus}}.tsv"
    params:
        configfile = CONFIGFILE
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} CHROMSIM_LOCUS={wildcards.locus} \
            "$CONDA_PREFIX/bin/python" scripts/validate_ctcf.py 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule simulateShard:
    input:
        validated = f"{L('done')}/validateCTCF.done",
        polysim   = rules.preparePolysim.output.done,
        bed       = lambda wildcards: rel(ctcf_bed_path(
            wildcards.locus, COND_LOOKUP[wildcards.condition]["ctcf_type"])),
        ## peak files for the condition's sticky / loading / pausing options (none by default)
        features  = lambda wildcards: [rel(f) for f in
                                       condition_peak_files(COND_LOOKUP[wildcards.condition]).values()]
    output:
        done = f"{L('done')}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.done"
    log:
        err = f"{LOG_DIR}/{{locus}}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.err",
        out = f"{LOG_DIR}/{{locus}}/simulateShard_{{condition}}_rep{{rep}}_shard{{shard}}.out"
    benchmark:
        f"{BENCH_DIR}/simulateShard_{{locus}}_{{condition}}_rep{{rep}}_shard{{shard}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = lambda wildcards: LP(wildcards, "sims"),
        n_shards   = N_SHARDS
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} CHROMSIM_LOCUS={wildcards.locus} \
            "$CONDA_PREFIX/bin/python" scripts/run_simulation_shard.py \
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
        lambda wildcards: expand(rules.simulateShard.output.done, locus=wildcards.locus,
                                 condition=wildcards.condition, rep=wildcards.rep,
                                 shard=range(N_SHARDS))
    output:
        done = f"{L('done')}/mergeShards_{{condition}}_rep{{rep}}.done"
    log:
        err = f"{LOG_DIR}/{{locus}}/mergeShards_{{condition}}_rep{{rep}}.err",
        out = f"{LOG_DIR}/{{locus}}/mergeShards_{{condition}}_rep{{rep}}.out"
    benchmark:
        f"{BENCH_DIR}/mergeShards_{{locus}}_{{condition}}_rep{{rep}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = lambda wildcards: LP(wildcards, "sims"),
        n_shards   = N_SHARDS
    conda: ENV_SIM
    shell:
        """
        CHROMSIM_CONFIG={params.configfile} CHROMSIM_LOCUS={wildcards.locus} \
            "$CONDA_PREFIX/bin/python" scripts/merge_shards.py \
            --condition {wildcards.condition} \
            --replicate {wildcards.rep} \
            --n-shards {params.n_shards} \
            --results-dir {params.sim_dir} 1> {log.out} 2> {log.err}
        touch {output.done}
        """


rule analyzeCondition:
    input:
        merged = lambda wildcards: expand(rules.mergeShards.output.done, locus=wildcards.locus,
                                          condition=wildcards.condition,
                                          rep=range(N_REPLICATES)),
        mcool  = compare_mcool,
        ## with calibration.msd.reference_condition, other conditions wait for its units
        reference = lambda wildcards: (
            [f"{OUT}/loci/{wildcards.locus}/done/analyzeCondition_{CALIBRATION['msd']['reference_condition']}.done"]
            if CALIBRATION["msd"]["reference_condition"]
            and CALIBRATION["default_anchor"] == "msd"
            and wildcards.condition != CALIBRATION["msd"]["reference_condition"] else [])
    output:
        done = f"{L('done')}/analyzeCondition_{{condition}}.done"
    log:
        err = f"{LOG_DIR}/{{locus}}/analyzeCondition_{{condition}}.err",
        out = f"{LOG_DIR}/{{locus}}/analyzeCondition_{{condition}}.out"
    benchmark:
        f"{BENCH_DIR}/analyzeCondition_{{locus}}_{{condition}}.tsv"
    params:
        configfile = CONFIGFILE,
        sim_dir    = lambda wildcards: LP(wildcards, "sims"),
        out_dir    = lambda wildcards: f"{LP(wildcards, 'analysis')}/{wildcards.condition}",
        mcool_dir  = MCOOL_DIR,
        n_jobs     = ANALYSIS["n_jobs"],
        msd_flag   = "" if ANALYSIS["msd"] else "--no-msd"
    conda: ENV_SIM
    shell:
        """
        mkdir -p {params.out_dir}
        CHROMSIM_CONFIG={params.configfile} CHROMSIM_LOCUS={wildcards.locus} \
            "$CONDA_PREFIX/bin/python" scripts/run_analysis_all.py \
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
        expand(rules.analyzeCondition.output.done, locus=LOCUS_NAMES, condition=COND_NAMES)
    output:
        f"{OUT}/{RUN_NAME}_benchmarking.tsv"
    log:
        err = f"{LOG_DIR}/benchmarkSummary.err"
    params:
        bench_dir = BENCH_DIR
    conda: ENV_SIM
    shell:
        """
        "$CONDA_PREFIX/bin/python" workflows/utils/benchmarking.py {params.bench_dir} {output} 2> {log.err}
        """
