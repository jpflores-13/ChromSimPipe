#!/bin/bash
#SBATCH --job-name=merge_shards
#SBATCH --output=logs/merge_shards_%j.out
#SBATCH --error=logs/merge_shards_%j.err
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
# Walltime pinned to the cluster's 48 h QOS ceiling. SLURM always
# requires a --time, so this is the closest equivalent to "no time
# limit" -- the job will run as long as it needs and only get killed
# if it actually exceeds 48 h. A single condition/replicate merge of 8
# polychrom shards with 25k frames each takes ~1 h; with 8 outer
# workers running condition groups in parallel, even a 30+ group fresh
# sweep finishes well inside this budget. (2 h budget was hit by job
# 9448832 on 2026-04-27 with TIMEOUT; 8 h was a safer interim; 48 h
# removes the worry entirely.)

# ============================================================================
# MERGE SHARDS for every condition / replicate present under
# results/polychrom_3d/. Produces directories named
#   results/polychrom_3d/merged_<condition>_rep<N>/
# (the 'merged_' prefix is the 2026-04-27 naming convention so merged
# raw-data dirs are immediately distinguishable from individual shard dirs).
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_merge_shards.sh
#
#   # remove shard dirs after a successful merge (saves disk):
#   CLEANUP=1 sbatch cluster/submit_merge_shards.sh
#
#   # disable parallelism (one group at a time, easier to debug):
#   NO_PARALLEL=1 sbatch cluster/submit_merge_shards.sh
#
# Pre-existing legacy un-prefixed merged dirs (e.g.
# 'mESC_baseline_rep0/' instead of 'merged_mESC_baseline_rep0/') from
# earlier runs are not touched. Run scripts/migrate_legacy_analysis.py
# with --rename-to-merged to bring them under the new naming.
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

ARGS=(--all --results-dir "${RESULTS_DIR}")

[ "${CLEANUP:-0}"     = "1" ] && ARGS+=(--cleanup)
[ "${NO_PARALLEL:-0}" = "1" ] && ARGS+=(--no-parallel)

if [ "$#" -gt 0 ]; then
    ARGS+=("$@")
fi

echo "============================================================"
echo "MERGE SHARDS - $(date)"
echo "============================================================"
echo "Job ID:     ${SLURM_JOB_ID:-<interactive>}"
echo "Node:       $(hostname)"
echo "Results:    ${RESULTS_DIR}"
echo "Cleanup:    ${CLEANUP:-0}"
echo "Args:       ${ARGS[*]}"
echo "============================================================"
echo ""

python scripts/merge_shards.py "${ARGS[@]}"

echo ""
echo "Done. Merged dirs are at ${RESULTS_DIR}/merged_*"
