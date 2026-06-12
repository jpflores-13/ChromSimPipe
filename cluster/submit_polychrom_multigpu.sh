#!/bin/bash
#SBATCH --job-name=poly3d
#SBATCH --output=logs/poly3d_%A_%a.out
#SBATCH --error=logs/poly3d_%A_%a.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:4
#SBATCH --time=48:00:00
#SBATCH --array=0-7       # 8 conditions, one job per condition

# ============================================================================
# MULTI-GPU POLYCHROM SIMULATION (4 V100s per condition)
#
# Strategy:
#   - Each SLURM array task = 1 condition (from SIMULATION_CONDITIONS)
#   - Within each task, we launch 4 sub-simulations in parallel (1 per GPU)
#   - Sub-simulations split the total blocks across GPUs
#   - After all finish, we merge conformations into a single trajectory
#
# This gives 4× speedup per condition within the 48-hour wall time.
#
# Total: 8 array jobs × 4 GPUs each = 32 GPU-hours per condition
#        (enough for ~20,000 blocks at 2000 monomers)
#
# Usage:
#   cd cohesin_sim && mkdir -p logs
#   sbatch cluster/submit_polychrom_multigpu.sh
#
# To run a subset:
#   sbatch --array=0,3 cluster/submit_polychrom_multigpu.sh  # mESC_ctrl + CN_long_residency
# ============================================================================

set -euo pipefail

# --- Load environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# --- Condition mapping ---
# These condition names match SIMULATION_CONDITIONS in parameters.py.
CONDITIONS=(
    "mESC_ctrl"                          # mESC params + mESC CTCF
    "mESC_params_neuron_ctcf"            # mESC params + neuron CTCF
    "CN_baseline_neuron_ctcf"            # null hypothesis
    "CN_long_residency_neuron_ctcf"      # 2× processivity
    "CN_very_long_residency_neuron_ctcf" # 4× processivity
    "CN_high_density_neuron_ctcf"        # 1.5× cohesin density
    "CN_long_res_high_dens_neuron_ctcf"  # 2× proc + 1.5× density
    "CN_long_res_low_dens_neuron_ctcf"   # 3× proc, lower density
)
CONDITION=${CONDITIONS[$SLURM_ARRAY_TASK_ID]}

N_GPUS=4
N_REPLICATES=3    # we run 3 full replicates, each split across GPUs
REP_OFFSET=3      # start at rep3 so we don't clash with single-GPU jobs (reps 0-2)

# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"

echo "============================================"
echo "Condition:     ${CONDITION}"
echo "GPUs:          ${N_GPUS}"
echo "Replicates:    ${N_REPLICATES}"
echo "Node:          $(hostname)"
echo "GPUs avail:    $(nvidia-smi -L 2>/dev/null | wc -l)"
echo "Date:          $(date)"
echo "============================================"

# --- Launch sub-simulations in parallel ---
# For each replicate, split total_blocks across N_GPUS shards.
# Each shard runs on a separate GPU and produces its own conformations.
# Then we merge all shards for that replicate.

RESULTS_BASE="results/polychrom_3d"
PIDS=()

for REP in $(seq ${REP_OFFSET} $((REP_OFFSET + N_REPLICATES - 1))); do
    for GPU_IDX in $(seq 0 $((N_GPUS - 1))); do
        SHARD_DIR="${RESULTS_BASE}/${CONDITION}_rep${REP}_shard${GPU_IDX}"

        echo "[$(date +%H:%M:%S)] Launching: ${CONDITION} rep${REP} shard${GPU_IDX} on GPU ${GPU_IDX}"

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
done

# Wait for all sub-simulations
echo "Waiting for ${#PIDS[@]} sub-simulations..."
FAIL=0
for PID in "${PIDS[@]}"; do
    wait "$PID" || FAIL=$((FAIL + 1))
done

if [ "$FAIL" -gt 0 ]; then
    echo "WARNING: ${FAIL} sub-simulations failed"
fi

# --- Merge shards ---
echo ""
echo "Merging shards..."
for REP in $(seq ${REP_OFFSET} $((REP_OFFSET + N_REPLICATES - 1))); do
    echo "[$(date +%H:%M:%S)] Merging: ${CONDITION} rep${REP}"
    python scripts/merge_shards.py \
        --condition "${CONDITION}" \
        --replicate "${REP}" \
        --n-shards "${N_GPUS}" \
        --results-dir "${RESULTS_BASE}"
done

echo ""
echo "All done: ${CONDITION} — $(date)"
