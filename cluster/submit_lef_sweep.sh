#!/bin/bash
#SBATCH --job-name=lef_sweep
#SBATCH --output=logs/lef_sweep_%A_%a.out
#SBATCH --error=logs/lef_sweep_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-23    # 8 conditions × 3 replicates = 24 jobs (0-23)

# ============================================================================
# LEF-ONLY PARAMETER SWEEP (CPU, no GPU needed)
#
# Runs all 8 simulation conditions × 3 replicates as a SLURM job array.
# Each condition pairs a cohesin parameter set with a CTCF site set
# (mESC or neuron), as defined in SIMULATION_CONDITIONS in parameters.py.
#
# Each job takes ~10-30 min on a single CPU core.
#
# Usage:
#   cd cohesin_sim
#   mkdir -p logs
#   sbatch cluster/submit_lef_sweep.sh
# ============================================================================

set -euo pipefail

# --- Load environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# --- Map array index to (condition, replicate) ---
# These condition names match SIMULATION_CONDITIONS in parameters.py.
# Each condition pairs a cohesin parameter set with the correct CTCF site set.
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
echo "Date:          $(date)"
echo "============================================"

# --- Run simulation ---
# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"

python scripts/run_simulation.py \
    --condition "${CONDITION}" \
    --replicate "${REPLICATE_IDX}" \
    --engine lef_only \
    --output results/lef_sweep

echo "Done: ${CONDITION} rep${REPLICATE_IDX} — $(date)"
