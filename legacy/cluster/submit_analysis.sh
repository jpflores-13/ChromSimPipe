#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH --output=logs/analysis_%A_%a.out
#SBATCH --error=logs/analysis_%A_%a.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-23    # 8 conditions × 3 replicates = 24 jobs (0-23)

# ============================================================================
# POST-SIMULATION ANALYSIS (single-dir, legacy)
#
# Computes contact map + P(s) + insulation for ONE condition/replicate at a
# time against experimental Hi-C. This is the simple legacy entry point that
# calls analysis/contact_maps.py directly on a single simulation directory.
#
# For the full pipeline (APA, CTCF overlays, mcool comparison, pooling
# across replicates, MSD/dynamics, calibration), prefer
#   cluster/submit_analysis_cpu.sh
# which calls scripts/run_analysis_all.py and walks ALL directories at once.
#
# Prerequisites:
#   - Simulation results in results/polychrom_3d/ (or results/lef_sweep/)
#   - Experimental Hi-C matrices in data/ (optional but recommended)
#
# Usage:
#   sbatch cluster/submit_analysis.sh
# ============================================================================

set -euo pipefail

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# Full condition names from SIMULATION_CONDITIONS in configs/parameters.py.
# These MUST match exactly — run_simulation.py writes output dirs as
#   <params_name>_ctcf-<ctcf_type>_rep<N>
# so we reconstruct that suffix from the condition metadata below.
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
N_REPLICATES=3

CONDITION_IDX=$(( SLURM_ARRAY_TASK_ID / N_REPLICATES ))
REPLICATE_IDX=$(( SLURM_ARRAY_TASK_ID % N_REPLICATES ))
CONDITION=${CONDITIONS[$CONDITION_IDX]}

# SLURM copies scripts to /var/spool, so $0 won't point to our repo.
cd "${SLURM_SUBMIT_DIR}"

# --- Choose which results to analyze ---
# Change this to results/lef_sweep if you ran LEF-only
RESULTS_DIR="results/polychrom_3d"

# Translate the condition name (as stored in SIMULATION_CONDITIONS) into the
# actual simulation output directory name. run_simulation.py builds the dir
# as "<params.name>_ctcf-<ctcf_type>_rep<N>"; we ask parameters.py itself to
# resolve (params.name, ctcf_type) so we stay in sync with any edits there.
read -r PARAMS_NAME CTCF_TYPE < <(python -c "
from configs.parameters import get_condition
c = get_condition('${CONDITION}')
print(c['params']['name'], c['ctcf_type'])
") || {
    echo "ERROR: could not resolve condition '${CONDITION}' via configs.parameters.get_condition"
    exit 1
}

SIM_DIR="${RESULTS_DIR}/${PARAMS_NAME}_ctcf-${CTCF_TYPE}_rep${REPLICATE_IDX}"

if [ ! -d "${SIM_DIR}" ]; then
    echo "No results found at ${SIM_DIR}, skipping."
    exit 0
fi

echo "Analyzing: ${CONDITION} rep${REPLICATE_IDX}  (dir: ${SIM_DIR})"

# --- Run contact map analysis ---
ARGS="--sim-dir ${SIM_DIR}"

# Add experimental Hi-C if available (edit paths as needed).
# The right reference depends on CTCF_TYPE: neuron sims compare to CN Hi-C.
if [ "${CTCF_TYPE}" = "neuron" ] && [ -f "data/hic_CN_Sox2.npy" ]; then
    ARGS="${ARGS} --exp-hic data/hic_CN_Sox2.npy"
elif [ -f "data/hic_mESC_Sox2.npy" ]; then
    ARGS="${ARGS} --exp-hic data/hic_mESC_Sox2.npy"
fi

python analysis/contact_maps.py ${ARGS}

echo "Done: ${CONDITION} rep${REPLICATE_IDX}"
