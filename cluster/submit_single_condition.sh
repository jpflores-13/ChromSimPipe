#!/bin/bash
#SBATCH --job-name=cohesin_1
#SBATCH --output=logs/cohesin_1_%j.out
#SBATCH --error=logs/cohesin_1_%j.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00

# ============================================================================
# RUN A SINGLE CONDITION ON 4× V100
#
# Useful for testing one condition before running all 8.
#
# Usage:
#   cd cohesin_sim && mkdir -p logs
#
#   # Run the mESC baseline:
#   sbatch cluster/submit_single_condition.sh mESC_ctrl
#
#   # Run a specific neuron hypothesis:
#   sbatch cluster/submit_single_condition.sh CN_long_residency_neuron_ctcf
#
# Available conditions (from parameters.py):
#   mESC_ctrl                          — mESC params + mESC CTCF
#   mESC_params_neuron_ctcf            — mESC params + neuron CTCF
#   CN_baseline_neuron_ctcf            — null hypothesis
#   CN_long_residency_neuron_ctcf      — 2× processivity
#   CN_very_long_residency_neuron_ctcf — 4× processivity
#   CN_high_density_neuron_ctcf        — 1.5× cohesin density
#   CN_long_res_high_dens_neuron_ctcf  — 2× proc + 1.5× density
#   CN_long_res_low_dens_neuron_ctcf   — 3× proc, lower density
# ============================================================================

set -euo pipefail

# --- Get condition from command line ---
# When submitting: sbatch cluster/submit_single_condition.sh <condition_name>
# SLURM passes extra arguments after the script name
CONDITION="${1:?Usage: sbatch cluster/submit_single_condition.sh <condition_name>}"

# --- Load environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

N_GPUS=4
N_REPLICATES=3
REP_OFFSET=3          # start at rep3 so we don't clash with single-GPU jobs (reps 0-2)
RESULTS_BASE="results/polychrom_3d"

# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"
mkdir -p "${RESULTS_BASE}" logs

echo "============================================================"
echo "SINGLE CONDITION: ${CONDITION}"
echo "Node:    $(hostname)"
echo "GPUs:    ${N_GPUS}× $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "Date:    $(date)"
echo "============================================================"

# Validate condition name
python -c "from configs.parameters import get_condition; get_condition('${CONDITION}')" || {
    echo "ERROR: Unknown condition '${CONDITION}'"
    echo "Available conditions:"
    python -c "from configs.parameters import list_conditions; list_conditions()"
    exit 1
}

for REP in $(seq ${REP_OFFSET} $((REP_OFFSET + N_REPLICATES - 1))); do
    echo ""
    echo "[$(date +%H:%M:%S)] Replicate ${REP} (of ${REP_OFFSET}-$((REP_OFFSET + N_REPLICATES - 1)))..."

    PIDS=()
    for GPU_IDX in $(seq 0 $((N_GPUS - 1))); do
        echo "  → shard ${GPU_IDX} on GPU ${GPU_IDX}"

        python scripts/run_simulation_shard.py \
            --condition "${CONDITION}" \
            --replicate "${REP}" \
            --shard-index "${GPU_IDX}" \
            --n-shards "${N_GPUS}" \
            --gpu "${GPU_IDX}" \
            --output "${RESULTS_BASE}" \
            &

        PIDS+=($!)
    done

    FAIL=0
    for PID in "${PIDS[@]}"; do
        wait "$PID" || FAIL=$((FAIL + 1))
    done

    if [ "$FAIL" -gt 0 ]; then
        echo "  WARNING: ${FAIL} shards failed"
    fi

    echo "[$(date +%H:%M:%S)] Merging shards..."
    python scripts/merge_shards.py \
        --condition "${CONDITION}" \
        --replicate "${REP}" \
        --n-shards "${N_GPUS}" \
        --results-dir "${RESULTS_BASE}"
done

echo ""
echo "Done: ${CONDITION} — $(date)"
echo "Results: ${RESULTS_BASE}/"
