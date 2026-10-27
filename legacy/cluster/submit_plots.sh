#!/bin/bash
#SBATCH --job-name=plots
#SBATCH --output=logs/plots_%j.out
#SBATCH --error=logs/plots_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# ============================================================================
# GENERATE FIGURES (single job, run after all analyses complete)
#
# Usage:
#   sbatch cluster/submit_plots.sh
#   # or with dependency on analysis jobs:
#   sbatch --dependency=afterok:<ANALYSIS_JOB_ID> cluster/submit_plots.sh
# ============================================================================

set -euo pipefail

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"

# Change to results/lef_sweep if you ran LEF-only
RESULTS_DIR="results/polychrom_3d"

ARGS="--results-dir ${RESULTS_DIR}"

if [ -d "data" ]; then
    ARGS="${ARGS} --data-dir data"
fi

python analysis/plot_results.py ${ARGS} --output "${RESULTS_DIR}/figures"

echo "Figures saved to ${RESULTS_DIR}/figures/"
