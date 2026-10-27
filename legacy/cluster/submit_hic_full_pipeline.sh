#!/bin/bash
#SBATCH --job-name=hic_full
#SBATCH --output=logs/hic_full_%j.out
#SBATCH --error=logs/hic_full_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# Bonev Hi-C: ONE sbatch to run the WHOLE pipeline (download -> extract).
#
# This is the "fire and forget" entrypoint. It runs all 4 steps serially
# from a single job, and self-resubmits if it approaches the 7-day walltime.
#
# Why this works for an absurdly long pipeline:
#   1. Every step is idempotent. If the job dies and restarts, each script
#      checks for its expected outputs and skips work that is already done:
#        - 00_download_fastq.sh checks per-SRR gz pairs and per-rep merged FASTQs
#        - 01_align_and_parse.sh checks for the dedup pairs.gz
#        - 02_make_mcool.sh checks for per-rep and merged mcools
#        - 03_extract_sox2.py overwrites .npy files (cheap)
#   2. SLURM is configured to send SIGUSR1 to this batch script 10 minutes
#      before the walltime is reached (--signal=B:USR1@600). The trap below
#      sets a flag; between major steps we check the flag and, if set,
#      resubmit a fresh copy of this script and exit cleanly. The new job
#      picks up exactly where this one left off (idempotent skip).
#   3. If a step crashes (non-zero exit), SLURM sends a FAIL email and the
#      job stops. There is no auto-retry on hard failure: the user must
#      diagnose. (Re-running the pipeline is then just resubmitting this
#      script; everything done so far is preserved.)
#
# What this script does NOT do:
#   - Run alignments in parallel. Each replicate is processed serially in
#     the loop below. If you have the cluster cycles to align all 8
#     replicates concurrently, prefer the array-based pipeline:
#         sbatch cluster/submit_hic_download.sh
#         sbatch --dependency=afterok:<DL>     cluster/submit_hic_pipeline.sh
#         sbatch --dependency=afterok:<ALIGN>  cluster/submit_hic_mcool.sh
#         sbatch --dependency=afterok:<MCOOL>  cluster/submit_hic_extract.sh
#     This orchestrator is for "I want to submit one thing and walk away".
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_hic_full_pipeline.sh
#
# Status:
#   tail -f logs/hic_full_*.out
#   ls -lh data/mcool/      # see which mcools have appeared
#   squeue -u $USER -n hic_full   # see if a successor was resubmitted
# ============================================================================

set -uo pipefail   # NOT -e: we handle errors ourselves so we can resubmit

THIS_SCRIPT="${SLURM_SUBMIT_DIR}/cluster/submit_hic_full_pipeline.sh"
SAMPLES=(ES_rep1 ES_rep2 ES_rep3 ES_rep4 CN_rep1 CN_rep2 CN_rep3 CN_rep4)

# Self-resubmit machinery: SIGUSR1 is sent by SLURM 10 minutes before
# walltime. Set a flag; we check between heavy steps.
TIMEOUT_NEAR=0
trap 'echo "[$(date +%H:%M:%S)] SIGUSR1 received, will resubmit after current step"; TIMEOUT_NEAR=1' USR1

resubmit_and_exit() {
    echo
    echo "[$(date)] Resubmitting orchestrator for the next chunk of work..."
    NEW_JOB=$(sbatch --parsable "${THIS_SCRIPT}")
    echo "    new job id: ${NEW_JOB}"
    echo "    exiting current orchestrator cleanly"
    exit 0
}

check_timeout() {
    if [ "${TIMEOUT_NEAR}" -eq 1 ]; then
        resubmit_and_exit
    fi
}

cd "${SLURM_SUBMIT_DIR}"

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-hic-bonev}"

echo "==============================================================="
echo " Bonev Hi-C full pipeline (single-sbatch orchestrator)"
echo " Job ID    : ${SLURM_JOB_ID}"
echo " Host      : $(hostname)"
echo " CPUs      : ${SLURM_CPUS_PER_TASK}"
echo " Walltime  : $(scontrol show job ${SLURM_JOB_ID} 2>/dev/null \
                       | awk -F= '/TimeLimit=/{print $4}' || echo unknown)"
echo " Start     : $(date)"
echo "==============================================================="

# === STEP 0: download ===
echo
echo "=== STEP 0: download FASTQs (8 replicates) ==="
if ! THREADS="${SLURM_CPUS_PER_TASK}" \
        bash scripts/process_hic/00_download_fastq.sh; then
    echo "ERROR: step 0 (download) failed. Stopping."
    exit 1
fi
check_timeout

# === STEP 1: align + parse + dedup, per replicate, idempotent ===
echo
echo "=== STEP 1: align + parse + dedup (per replicate) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    OUT="data/pairs/${SAMPLE}.dedup.pairs.gz"
    if [ -f "${OUT}" ]; then
        echo "  ${SAMPLE}: already done (${OUT}), skipping"
        continue
    fi

    echo
    echo "  --- ${SAMPLE} starting at $(date) ---"
    if ! NTHREADS="${SLURM_CPUS_PER_TASK}" \
            bash scripts/process_hic/01_align_and_parse.sh "${SAMPLE}"; then
        echo "ERROR: step 1 failed for ${SAMPLE}. Stopping."
        exit 1
    fi
    echo "  --- ${SAMPLE} done at $(date) ---"

    check_timeout
done

# === STEP 2: merge + mcool ===
echo
echo "=== STEP 2: merge + mcool ==="
if [ -f "data/mcool/ES.mcool" ] && [ -f "data/mcool/CN.mcool" ]; then
    echo "  merged mcools already present, skipping step 2"
else
    if ! NTHREADS="${SLURM_CPUS_PER_TASK}" \
            bash scripts/process_hic/02_make_mcool.sh; then
        echo "ERROR: step 2 (mcool) failed. Stopping."
        exit 1
    fi
fi
check_timeout

# === STEP 3: extract Sox2 locus ===
echo
echo "=== STEP 3: extract Sox2 ==="
if ! python scripts/process_hic/03_extract_sox2.py \
        --mcool-es data/mcool/ES.mcool \
        --mcool-cn data/mcool/CN.mcool \
        --output data; then
    echo "ERROR: step 3 (Sox2 extraction) failed. Stopping."
    exit 1
fi

echo
echo "==============================================================="
echo " PIPELINE COMPLETE at $(date)"
echo " Outputs:"
echo "   data/mcool/ES.mcool"
echo "   data/mcool/CN.mcool"
echo "   data/mcool/per_rep/*.mcool   (8 per-replicate mcools)"
echo "   data/hic_*.npy               (Sox2 matrices, P(s), insulation, ...)"
echo "==============================================================="
touch "${SLURM_SUBMIT_DIR}/.hic_pipeline_complete"
