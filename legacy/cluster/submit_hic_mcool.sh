#!/bin/bash
#SBATCH --job-name=hic_mcool
#SBATCH --output=logs/hic_mcool_%j.out
#SBATCH --error=logs/hic_mcool_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# STEP 2: Build per-replicate AND merged-per-cell-type mcool files
#
# Single job that runs scripts/process_hic/02_make_mcool.sh, which:
#   1. Builds and balances one mcool per replicate (data/mcool/per_rep/*.mcool)
#   2. Merges per-replicate dedup pairs into per-cell-type pairs
#   3. Builds and balances one merged mcool per cell type (data/mcool/*.mcool)
#
# Idempotent: skips any output that already exists.
#
# Walltime sized for the deep Bonev libraries; pairs-level merge of 4 reps
# can take 2-4 hours and ICE balancing across 10 resolutions takes another
# 4-8 hours. 24h gives headroom; bump if you adopt finer base resolution.
#
# Usage:
#   ALIGN_JOB=$(sbatch --parsable cluster/submit_hic_pipeline.sh)
#   sbatch --dependency=afterok:${ALIGN_JOB} cluster/submit_hic_mcool.sh
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-hic-bonev}"

echo "============================================"
echo " Bonev Hi-C mcool generation"
echo " Job ID      : ${SLURM_JOB_ID}"
echo " Node        : $(hostname)"
echo " CPUs        : ${SLURM_CPUS_PER_TASK}"
echo " Start       : $(date)"
echo "============================================"

NTHREADS="${SLURM_CPUS_PER_TASK}" \
    bash scripts/process_hic/02_make_mcool.sh

echo
echo "mcool generation complete at $(date)"
