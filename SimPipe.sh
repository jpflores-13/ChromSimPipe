#!/usr/bin/env bash
# ==============================================================================
# filename:    SimPipe.sh
# project:     STRSsim
# description: sbatch entry point for the STRSsim Snakemake pipeline.
#              Submits Snakemake as a SLURM job; Snakemake then launches all
#              downstream rules via the SLURM executor plugin.
#
# Usage:
#   sbatch SimPipe.sh
#
# Requirements (run setup_data.sh first):
#   - conda envs: cohesin_sim, ctcf_extraction
#   - data/mcool/STRS_control.mcool  (or let the pipeline convert it)
#   - data/ctcf_beds/*.bed            (or let the pipeline extract them)
#   - config/SimConfig.yaml           (edit paths before running)
#
# To run a dry-run first (shows all jobs without submitting):
#   bash SimPipe.sh --dry-run
# ==============================================================================

#SBATCH --job-name=STRSsim
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --partition=general

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

mkdir -p logs logs_slurm

# ── Parse optional --dry-run flag ─────────────────────────────────────────────
DRY_RUN=""
for arg in "$@"; do
    [[ "${arg}" == "--dry-run" ]] && DRY_RUN="--dry-run --printshellcmds"
done

# ── Set up Snakemake Python venv ──────────────────────────────────────────────
# We use a local venv for Snakemake itself (not a conda env) to match the
# bagPipes RHEL9 pattern: module load python + pip install snakemake.
# The simulation and analysis envs (cohesin_sim, ctcf_extraction) are conda
# envs invoked via `conda run` inside each rule's shell block.

PYTHON_MODULE="python/${PYTHONVERS:-3.12.4}"
VENV_DIR="${REPO_ROOT}/.snakemake_venv"

module load "${PYTHON_MODULE}" 2>/dev/null || module load python/3.12.4

if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
    echo "Creating Snakemake venv at ${VENV_DIR} ..."
    python -m venv "${VENV_DIR}"
    source "${VENV_DIR}/bin/activate"
    pip install --quiet --upgrade pip
    pip install --quiet \
        "snakemake==8.27.1" \
        "snakemake-executor-plugin-slurm" \
        "pandas"
    echo "Venv created."
else
    source "${VENV_DIR}/bin/activate"
fi

echo "Snakemake version: $(snakemake --version)"
echo "Active locus: $(grep '^ACTIVE_LOCUS' configs/parameters.py)"
echo ""

# ── Launch pipeline ────────────────────────────────────────────────────────────
snakemake \
    --snakefile    workflows/STRSsim.snakefile \
    --profile      profiles/slurm \
    --configfile   config/SimConfig.yaml \
    --jobs         100 \
    --rerun-incomplete \
    --latency-wait 500 \
    --keep-going \
    ${DRY_RUN}

echo ""
echo "Pipeline complete."
