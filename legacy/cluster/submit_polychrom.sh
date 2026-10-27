#!/bin/bash
#SBATCH --job-name=polychrom
#SBATCH --output=logs/polychrom_%A_%a.out
#SBATCH --error=logs/polychrom_%A_%a.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --array=0-23    # 8 conditions × 3 replicates = 24 jobs (0-23)

# ============================================================================
# FULL 3D POLYCHROM SIMULATION (GPU, V100)
#
# Runs loop extrusion + 3D polymer dynamics on GPU.
# Each job takes ~6-12 hours on a single V100 with default block count (5000).
#
# Each condition pairs a cohesin parameter set with the correct CTCF site set
# (mESC or neuron), as defined in SIMULATION_CONDITIONS in parameters.py.
#
# Usage:
#   cd cohesin_sim
#   mkdir -p logs
#   sbatch cluster/submit_polychrom.sh
#
# To run a subset (e.g., only mESC_ctrl + CN_long_residency, 3 reps each):
#   sbatch --array=0-2,9-11 cluster/submit_polychrom.sh
# ============================================================================

set -euo pipefail

# --- Load environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# --- Map array index to (condition, replicate) ---
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
N_CONDITIONS=${#CONDITIONS[@]}
N_REPLICATES=3

CONDITION_IDX=$(( SLURM_ARRAY_TASK_ID / N_REPLICATES ))
REPLICATE_IDX=$(( SLURM_ARRAY_TASK_ID % N_REPLICATES ))
CONDITION=${CONDITIONS[$CONDITION_IDX]}

echo "============================================"
echo "Job array ID:  ${SLURM_ARRAY_TASK_ID}"
echo "Condition:     ${CONDITION}"
echo "Replicate:     ${REPLICATE_IDX}"
echo "Node:          $(hostname)"
echo "GPU:           $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "CUDA:          $(nvcc --version 2>/dev/null | tail -1 || echo 'N/A')"
echo "Date:          $(date)"
echo "============================================"

# --- Run simulation ---
# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
# SLURM_SUBMIT_DIR is the directory where sbatch was called.
cd "${SLURM_SUBMIT_DIR}"

python scripts/run_simulation.py \
    --condition "${CONDITION}" \
    --replicate "${REPLICATE_IDX}" \
    --engine polychrom \
    --gpu 0 \
    --output results/polychrom_3d

echo "Done: ${CONDITION} rep${REPLICATE_IDX} — $(date)"
