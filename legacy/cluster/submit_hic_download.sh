#!/bin/bash
#SBATCH --job-name=hic_dl
#SBATCH --output=logs/hic_dl_%A_%a.out
#SBATCH --error=logs/hic_dl_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=0-7
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# STEP 0: Download Bonev 2017 Hi-C FASTQ files from SRA (8-element array)
#
# One array task per biological replicate (ES_rep1..4, CN_rep1..4).
# Each task runs prefetch + fasterq-dump over every SRR for its assigned
# replicate, then concatenates them into a single (R1, R2) FASTQ pair.
#
# Disk footprint per replicate is roughly 100-500 GB of compressed FASTQ;
# total across all 8 is 3-5 TB.
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_hic_download.sh                       # all 8
#   sbatch --array=0-3 cluster/submit_hic_download.sh           # only ES
#   sbatch --array=0-7%2 cluster/submit_hic_download.sh         # 2 at a time
#
# Caveats:
#   - SRA throughput is bandwidth-bound. Running all 8 in parallel may
#     saturate the cluster's egress. Use %2 to throttle.
#   - SRA can be flaky. The 00_download_fastq.sh script is idempotent at
#     the per-run level, so resubmitting a failed task just resumes.
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
echo " Bonev Hi-C download"
echo " Sample      : ${SAMPLE}"
echo " Array idx   : ${SLURM_ARRAY_TASK_ID}"
echo " Job ID      : ${SLURM_JOB_ID}"
echo " Node        : $(hostname)"
echo " CPUs        : ${SLURM_CPUS_PER_TASK}"
echo " Start       : $(date)"
echo "============================================"

THREADS="${SLURM_CPUS_PER_TASK}" \
    bash scripts/process_hic/00_download_fastq.sh "${SAMPLE}"

echo
echo "Download complete: ${SAMPLE} at $(date)"
