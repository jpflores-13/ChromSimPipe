#!/bin/bash
#SBATCH --job-name=cohesin_sim
#SBATCH --output=logs/cohesin_sim_%j.out
#SBATCH --error=logs/cohesin_sim_%j.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --gres=gpu:2
#SBATCH --time=48:00:00

# ============================================================================
# COHESIN LOOP EXTRUSION SIMULATION — ALL CONDITIONS ON 4× V100
#
# Runs all 8 simulation conditions sequentially, using all 4 GPUs in parallel
# for each condition. Each condition runs 3 replicates, and within each
# replicate the simulation blocks are split across the 4 GPUs (shards).
#
# SAFE TO RE-RUN: each execution auto-detects existing shard directories and
# assigns the next available shard indices, so no previous data is ever
# overwritten. Re-running adds more independent conformations that are all
# pooled together at the merge step.
#
# First run:  reps 3-5, shard0 + shard1  → ~2 000 frames/replicate
# Second run: reps 3-5, shard2 + shard3  → ~4 000 frames/replicate total
# Third run:  reps 3-5, shard4 + shard5  → ~6 000 frames/replicate total
# (and so on)
#
# Usage:
#   cd cohesin_sim
#   mkdir -p logs
#   sbatch cluster/submit_all_4gpu.sh
#
# To run only specific conditions, edit the CONDITIONS array below or use:
#   sbatch cluster/submit_single_condition.sh <condition_name>
# ============================================================================

set -euo pipefail

# --- Load environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# --- Configuration ---
N_GPUS=2
N_REPLICATES=3
REP_OFFSET=3          # start at rep3 so we don't clash with single-GPU jobs (reps 0-2)
RESULTS_BASE="results/polychrom_3d"

# All 8 conditions from SIMULATION_CONDITIONS in parameters.py.
# Comment out any conditions you don't want to run.
CONDITIONS=(
    "mESC_ctrl"                          # mESC params + mESC CTCF (baseline)
    "mESC_params_neuron_ctcf"            # mESC params + neuron CTCF (CTCF-only effect)
    "CN_baseline_neuron_ctcf"            # null hypothesis (same cohesin, neuron CTCF)
    "CN_long_residency_neuron_ctcf"      # 2× processivity
    "CN_very_long_residency_neuron_ctcf" # 4× processivity
    "CN_high_density_neuron_ctcf"        # 1.5× cohesin density
    "CN_long_res_high_dens_neuron_ctcf"  # 2× processivity + 1.5× density
    "CN_long_res_low_dens_neuron_ctcf"   # 3× processivity, lower density
)

# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"
mkdir -p "${RESULTS_BASE}" logs

# ---------------------------------------------------------------------------
# Helper: find the next available shard index for a condition + replicate.
#
# Looks for existing shard directories matching the condition/replicate prefix
# and returns (highest existing shard index + 1).  Returns 0 if none exist.
#
# Usage: OFFSET=$(get_shard_offset "mESC_ctrl" 3)
# ---------------------------------------------------------------------------
get_shard_offset() {
    local condition="$1"
    local rep="$2"
    python - <<EOF
import os, sys
sys.path.insert(0, '.')
try:
    from configs.parameters import get_condition
    cond   = get_condition('${condition}')
    prefix = "{}_ctcf-{}_rep${rep}".format(cond['params']['name'], cond['ctcf_type'])
    results = '${RESULTS_BASE}'
    max_shard = -1
    if os.path.isdir(results):
        for entry in os.listdir(results):
            if entry.startswith(prefix + '_shard'):
                try:
                    idx = int(entry.split('_shard')[-1])
                    max_shard = max(max_shard, idx)
                except ValueError:
                    pass
    print(max_shard + 1)
except Exception as e:
    print(0, file=sys.stderr)
    print(0)
EOF
}

# --- Print job info ---
echo "============================================================"
echo "COHESIN SIMULATION — $(date)"
echo "============================================================"
echo "Node:          $(hostname)"
echo "GPUs:          ${N_GPUS}× $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "CUDA:          $(nvcc --version 2>/dev/null | grep release | awk '{print $6}' || echo 'N/A')"
echo "Conditions:    ${#CONDITIONS[@]}"
echo "Replicates:    ${N_REPLICATES} per condition (reps ${REP_OFFSET}-$((REP_OFFSET + N_REPLICATES - 1)))"
echo "Results:       ${RESULTS_BASE}"
echo "NOTE: Existing shards will be preserved; new shards will be added."
echo "============================================================"
echo ""

# --- List conditions to run ---
python -c "
from configs.parameters import list_conditions
list_conditions()
"
echo ""

# --- Run each condition ---
TOTAL=${#CONDITIONS[@]}
FAILED_CONDITIONS=()

for COND_IDX in $(seq 0 $((TOTAL - 1))); do
    CONDITION=${CONDITIONS[$COND_IDX]}

    echo ""
    echo "============================================================"
    echo "[$(date +%H:%M:%S)] CONDITION $((COND_IDX + 1))/${TOTAL}: ${CONDITION}"
    echo "============================================================"

    for REP in $(seq ${REP_OFFSET} $((REP_OFFSET + N_REPLICATES - 1))); do

        # Detect the next available shard index so we never overwrite existing data
        SHARD_OFFSET=$(get_shard_offset "${CONDITION}" "${REP}")
        echo "[$(date +%H:%M:%S)] rep${REP}: starting at shard index ${SHARD_OFFSET} (${N_GPUS} new shards)"

        PIDS=()
        for GPU_IDX in $(seq 0 $((N_GPUS - 1))); do
            GLOBAL_SHARD=$((SHARD_OFFSET + GPU_IDX))
            echo "  → shard ${GLOBAL_SHARD} on GPU ${GPU_IDX}"

            python scripts/run_simulation_shard.py \
                --condition "${CONDITION}" \
                --replicate "${REP}" \
                --shard-index "${GPU_IDX}" \
                --shard-index-offset "${SHARD_OFFSET}" \
                --n-shards "${N_GPUS}" \
                --gpu "${GPU_IDX}" \
                --output "${RESULTS_BASE}" \
                &

            PIDS+=($!)
        done

        # Wait for all shards to finish
        SHARD_FAIL=0
        for PID in "${PIDS[@]}"; do
            wait "$PID" || SHARD_FAIL=$((SHARD_FAIL + 1))
        done

        if [ "$SHARD_FAIL" -gt 0 ]; then
            echo "  WARNING: ${SHARD_FAIL}/${N_GPUS} shards failed for ${CONDITION} rep${REP}"
            FAILED_CONDITIONS+=("${CONDITION}_rep${REP}")
        fi

        # Merge all shards found for this replicate (--all auto-discovers new ones)
        echo "[$(date +%H:%M:%S)] Merging all shards for ${CONDITION} rep${REP}..."
        python scripts/merge_shards.py \
            --all \
            --results-dir "${RESULTS_BASE}" \
            --no-parallel   # one condition at a time here; no need for intra-merge parallelism

    done

    echo "[$(date +%H:%M:%S)] Done: ${CONDITION}"
done

# --- Summary ---
echo ""
echo "============================================================"
echo "ALL DONE — $(date)"
echo "============================================================"
echo "Completed: ${TOTAL} conditions × ${N_REPLICATES} replicates"
echo "Results:   ${RESULTS_BASE}/"

if [ ${#FAILED_CONDITIONS[@]} -gt 0 ]; then
    echo ""
    echo "WARNING — The following condition/replicates had shard failures:"
    for F in "${FAILED_CONDITIONS[@]}"; do
        echo "  - ${F}"
    done
fi

echo ""
echo "Next steps:"
echo "  1. Check shard counts:  ls ${RESULTS_BASE}/ | grep shard | sort"
echo "  2. Run analysis:        python scripts/run_analysis_all.py \\"
echo "         --results-dir ${RESULTS_BASE} --skip-existing"
echo "  3. (Optional) Free disk space after verifying analysis:"
echo "         python scripts/merge_shards.py --all \\"
echo "             --results-dir ${RESULTS_BASE} --cleanup"
