#!/bin/bash
#SBATCH --job-name=migrate_legacy
#SBATCH --output=logs/migrate_legacy_%j.out
#SBATCH --error=logs/migrate_legacy_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00

# ============================================================================
# MIGRATE pre-2026-04-27 per-sim-dir analysis outputs to the flat
# results/analysis/ layout.
#
# What this does
# --------------
# For every results/polychrom_3d/<sim_dir>/ that has an analysis/
# subdirectory (the old layout), this:
#   1. moves each file to results/analysis/{condition}_*blk_rep*_*
#      using the same renaming rule the live pipeline applies;
#   2. removes the now-empty <sim_dir>/analysis/ directory;
#   3. (optional) renames <sim_dir> to merged_<sim_dir> so the new
#      naming convention applies.
#
# Usage
# -----
#   mkdir -p logs
#
#   # 1. dry-run first (recommended): logs intent without touching files
#   DRY_RUN=1 sbatch cluster/submit_migrate_legacy.sh
#
#   # 2. live run, move files, delete empty legacy dirs, rename to merged_
#   sbatch cluster/submit_migrate_legacy.sh
#
#   # 3. variants
#   COPY=1 sbatch cluster/submit_migrate_legacy.sh             # copy instead of move
#   NO_RENAME=1 sbatch cluster/submit_migrate_legacy.sh        # keep old dir names
#   NO_DELETE_EMPTY=1 sbatch cluster/submit_migrate_legacy.sh  # keep empty analysis/ dirs
#
# Defaults: live move, --delete-empty, --rename-to-merged.
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

RESULTS_DIR="results/polychrom_3d"
ANALYSIS_DIR="results/analysis"

ARGS=(--results-dir "${RESULTS_DIR}" --analysis-dir "${ANALYSIS_DIR}")

[ "${DRY_RUN:-0}"         = "1" ] && ARGS+=(--dry-run)
[ "${COPY:-0}"            = "1" ] && ARGS+=(--copy)
[ "${NO_DELETE_EMPTY:-0}" != "1" ] && ARGS+=(--delete-empty)
[ "${NO_RENAME:-0}"       != "1" ] && ARGS+=(--rename-to-merged)

if [ "$#" -gt 0 ]; then
    ARGS+=("$@")
fi

echo "============================================================"
echo "MIGRATE LEGACY ANALYSIS - $(date)"
echo "============================================================"
echo "Job ID:     ${SLURM_JOB_ID:-<interactive>}"
echo "Source:     ${RESULTS_DIR}"
echo "Dest:       ${ANALYSIS_DIR}"
echo "Mode:       $([ "${DRY_RUN:-0}" = "1" ] && echo dry-run || ([ "${COPY:-0}" = "1" ] && echo copy || echo move))"
echo "Delete empty <sim_dir>/analysis/: ${NO_DELETE_EMPTY:-0} (0 = yes)"
echo "Rename to merged_:                ${NO_RENAME:-0} (0 = yes)"
echo "Args:       ${ARGS[*]}"
echo "============================================================"
echo ""

python scripts/migrate_legacy_analysis.py "${ARGS[@]}"

echo ""
echo "Done. Files now under ${ANALYSIS_DIR}/"
