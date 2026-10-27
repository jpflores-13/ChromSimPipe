#!/bin/bash
# ============================================================================
# SLURM wrapper for scripts/download_bonev_hic.sh
# ----------------------------------------------------------------------------
# Submit ES and CN as two independent jobs so they run in parallel.
#
# Usage:
#   mkdir -p logs
#   sbatch --export=CELLTYPE=ES scripts/sbatch_bonev_hic.sbatch
#   sbatch --export=CELLTYPE=CN scripts/sbatch_bonev_hic.sbatch
#
# Or as a job array (ES = task 0, CN = task 1):
#   sbatch --array=0-1 scripts/sbatch_bonev_hic.sbatch
#
# Notes:
#   • This is CPU work (cooler cload + balance), not GPU.  Pick a CPU
#     partition — the `cuda` partition works but wastes a V100 slot.
#   • --cpus-per-task propagates into the script as NPROC via
#     $SLURM_CPUS_PER_TASK, driving both --nproc and --max-split.
#     With 40 cores chr3 is split into 40 chunks and every worker stays
#     busy during binning — expect an ~8-10× speed-up over the 4-core
#     local default.
#   • Memory: chr3 at 1 kb × ~200 M pairs peaks around 12 GB for cload
#     and ~1-2 GB per worker for balance.  Bumping to 40 workers means
#     budgeting ~80 GB peak, so we ask for 96 GB to have headroom.
# ============================================================================

#SBATCH --job-name=bonev_hic
#SBATCH --output=logs/bonev_%x_%A_%a.out
#SBATCH --error=logs/bonev_%x_%A_%a.err
#SBATCH --cpus-per-task=40
#SBATCH --mem=96G
#SBATCH --time=03:00:00
# Uncomment and set to your CPU partition if you have one:
# #SBATCH --partition=workq

set -euo pipefail

# --- Resolve which cell type to run --------------------------------------
if [ -n "${CELLTYPE:-}" ]; then
    CT="${CELLTYPE}"
elif [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CTS=(ES CN)
    CT="${CTS[$SLURM_ARRAY_TASK_ID]}"
else
    echo "ERROR: set CELLTYPE (--export=CELLTYPE=ES) or use --array=0-1" >&2
    exit 1
fi

echo "================================================================"
echo " Bonev Hi-C download — cell type: ${CT}"
echo " Host       : $(hostname)"
echo " CPUs       : ${SLURM_CPUS_PER_TASK:-?}"
echo " Submit dir : ${SLURM_SUBMIT_DIR:-$(pwd)}"
echo " Job ID     : ${SLURM_JOB_ID:-local}"
echo "================================================================"

# --- Activate the project conda env --------------------------------------
# No `module` system on this cluster — we add mambaforge to PATH directly,
# then use `eval "$(conda shell.bash hook)" + conda activate`.  Do NOT
# use `source activate`: on compute nodes it updates the prompt but
# silently keeps base Python 3.10 instead of the env's Python 3.12.
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# --- Run -----------------------------------------------------------------
bash scripts/download_bonev_hic.sh "${CT}"
