#!/bin/bash
#SBATCH --job-name=hic_extract
#SBATCH --output=logs/hic_extract_%j.out
#SBATCH --error=logs/hic_extract_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL

# ============================================================================
# STEP 3: Extract Sox2 locus matrices and derived quantities from mcools
#
# Runs scripts/process_hic/03_extract_sox2.py to produce per-cell-type:
#   - balanced and raw contact matrices at the Sox2 region
#   - O/E and expected matrices
#   - P(s) curves
#   - insulation profiles at multiple windows
#   - compartment eigenvector at 25 kb
#
# Idempotent: writes .npy files; rerun overwrites them.
#
# Usage:
#   MCOOL_JOB=$(sbatch --parsable cluster/submit_hic_mcool.sh)
#   sbatch --dependency=afterok:${MCOOL_JOB} cluster/submit_hic_extract.sh
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV:-hic-bonev}"

echo "============================================"
echo " Bonev Hi-C: Sox2 locus extraction"
echo " Job ID      : ${SLURM_JOB_ID}"
echo " Node        : $(hostname)"
echo " CPUs        : ${SLURM_CPUS_PER_TASK}"
echo " Start       : $(date)"
echo "============================================"

python scripts/process_hic/03_extract_sox2.py \
    --mcool-es data/mcool/ES.mcool \
    --mcool-cn data/mcool/CN.mcool \
    --output data

echo
echo "Sox2 extraction complete at $(date)"
