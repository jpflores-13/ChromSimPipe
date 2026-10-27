#!/bin/bash
#SBATCH --job-name=comparable_merge
#SBATCH --output=logs/comparable_merge_%j.out
#SBATCH --error=logs/comparable_merge_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# submit_comparable_merge.sh — materialise downsampled rep dirs for
# cross-condition comparable analysis (200 000-frame cap per condition).
#
# Wraps scripts/build_comparable_merge.py. Defaults to --reps 3
# --blocks-per-rep 700 --mode tail. Override via env vars or by passing
# extra args after the script name.
#
# Also clears stale `*_pooled_*` analysis artifacts so the next analysis
# run rebuilds them from current rep state (otherwise SKIP_EXISTING
# inside the analysis pipeline would short-circuit — see CLAUDE.md
# "Pooled analysis artifacts go stale silently").
#
# Usage
# -----
#   # As an auto-chained step from submit_catch_up.sh planner — handled
#   # automatically when you run `bash cluster/submit_catch_up.sh`.
#
#   # Standalone (e.g. after a manual catch-up round):
#   sbatch cluster/submit_comparable_merge.sh
#
#   # Override the comparable-set parameters:
#   COMPARABLE_REPS=5 COMPARABLE_BLOCKS=500 sbatch cluster/submit_comparable_merge.sh
#
#   # Skip the stale-pooled cleanup (e.g. if you only want comparable dirs
#   # without invalidating existing pooled artifacts):
#   CLEAN_STALE_POOLS=0 sbatch cluster/submit_comparable_merge.sh
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

COMPARABLE_REPS="${COMPARABLE_REPS:-3}"
COMPARABLE_BLOCKS="${COMPARABLE_BLOCKS:-700}"
COMPARABLE_MODE="${COMPARABLE_MODE:-tail}"
CLEAN_STALE_POOLS="${CLEAN_STALE_POOLS:-1}"

echo "============================================================"
echo "comparable_merge step  ($(date))"
echo "  reps/cond:    ${COMPARABLE_REPS}"
echo "  blocks/rep:   ${COMPARABLE_BLOCKS}"
echo "  mode:         ${COMPARABLE_MODE}"
echo "  clean pools:  ${CLEAN_STALE_POOLS}"
echo "============================================================"

python scripts/build_comparable_merge.py \
    --reps "${COMPARABLE_REPS}" \
    --blocks-per-rep "${COMPARABLE_BLOCKS}" \
    --mode "${COMPARABLE_MODE}"

if [ "${CLEAN_STALE_POOLS}" = "1" ]; then
    echo
    echo "Removing stale *_pooled_* artifacts so the next analysis rebuilds"
    echo "them from current rep state..."
    rm -fv results/analysis/*_pooled_* 2>/dev/null || true
    rm -fv results/analysis/*_rep*_contact_map.npy 2>/dev/null || true
fi

echo
echo "comparable_merge done — $(date)"
