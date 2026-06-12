#!/bin/bash
#SBATCH --job-name=plot_panels
#SBATCH --output=logs/plot_panels_%j.out
#SBATCH --error=logs/plot_panels_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# ============================================================================
# REGENERATE CONTACT-MAP FIGURE PANELS from results/analysis/.
#
# Produces three panel families under results/figures/:
#   * <condition>_contact_maps_per_rep.png       - one subplot per replicate
#   * <condition>_contact_map_merged_reps.png    - single map from pooled data
#   * all_conditions_contact_maps.png            - one subplot per condition
#
# This is what scripts/plot_contact_map_panels.py emits. The 16-CPU and
# 32-CPU analysis sbatch scripts already invoke this internally at the end
# of a run; use this standalone wrapper to regenerate panels without
# re-running the full analysis pipeline.
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_plot_panels.sh
#
#   # raw scale instead of log:
#   NO_LOG_SCALE=1 sbatch cluster/submit_plot_panels.sh
#
#   # chain after an analysis job:
#   sbatch --dependency=afterok:<ANALYSIS_JOB_ID> cluster/submit_plot_panels.sh
# ============================================================================

set -euo pipefail

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

if [ "${CONDA_DEFAULT_ENV:-}" != "polychrom" ]; then
    echo "ERROR: conda env is '${CONDA_DEFAULT_ENV:-<none>}', expected 'polychrom'." >&2
    exit 1
fi

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p logs

ANALYSIS_DIR="results/analysis"
OUTPUT_DIR="results/figures"

ARGS=(--analysis-dir "${ANALYSIS_DIR}" --output-dir "${OUTPUT_DIR}")
[ "${NO_LOG_SCALE:-0}" = "1" ] && ARGS+=(--no-log-scale)

if [ "$#" -gt 0 ]; then
    ARGS+=("$@")
fi

echo "============================================================"
echo "PLOT CONTACT-MAP PANELS - $(date)"
echo "============================================================"
echo "Analysis:   ${ANALYSIS_DIR}"
echo "Figures:    ${OUTPUT_DIR}"
echo "Args:       ${ARGS[*]}"
echo "============================================================"
echo ""

python scripts/plot_contact_map_panels.py "${ARGS[@]}"

echo ""
echo "Done. Panels in ${OUTPUT_DIR}/"
