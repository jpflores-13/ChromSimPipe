#!/bin/bash
#SBATCH --job-name=hic_align
#SBATCH --output=logs/hic_align_%A_%a.out
#SBATCH --error=logs/hic_align_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-7
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# STEP 1: Align + parse + dedup, one biological replicate per array task
#
# 8 array tasks: ES_rep1..4 + CN_rep1..4. Each runs the full
# bwa-mem2 mem -SP | pairtools parse | sort | dedup pipeline for one
# replicate at a time. Per-replicate dedup pairs land in data/pairs/.
#
# Idempotent: if data/pairs/<sample>.dedup.pairs.gz already exists, the
# task exits early. So resubmitting the same array on a partially-done
# dataset only redoes the missing replicates.
#
# Usage:
#   mkdir -p logs
#
#   # download first if not already done
#   sbatch cluster/submit_hic_download.sh
#
#   # then align (this script). Capture jobid for chaining.
#   ALIGN_JOB=$(sbatch --parsable cluster/submit_hic_pipeline.sh)
#
#   # then mcool (waits for all alignments to finish)
#   MCOOL_JOB=$(sbatch --parsable --dependency=afterok:${ALIGN_JOB} \
#                       cluster/submit_hic_mcool.sh)
#
#   # then extract
#   sbatch --dependency=afterok:${MCOOL_JOB} cluster/submit_hic_extract.sh
# ============================================================================

set -euo pipefail

SAMPLES=(
    "ES_rep1" "ES_rep2" "ES_rep3" "ES_rep4"
    "CN_rep1" "CN_rep2" "CN_rep3" "CN_rep4"
)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

cd "${SLURM_SUBMIT_DIR}"

# Activate conda env (override default with CONDA_ENV=<name> on sbatch)
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-hic-bonev}"

echo "============================================"
echo " Bonev Hi-C alignment + parsing"
echo " Sample      : ${SAMPLE}"
echo " Array idx   : ${SLURM_ARRAY_TASK_ID}"
echo " Job ID      : ${SLURM_JOB_ID}"
echo " Node        : $(hostname)"
echo " CPUs        : ${SLURM_CPUS_PER_TASK}"
echo " Start       : $(date)"
echo "============================================"

NTHREADS="${SLURM_CPUS_PER_TASK}" \
    bash scripts/process_hic/01_align_and_parse.sh "${SAMPLE}"

echo
echo "Alignment complete: ${SAMPLE} at $(date)"
