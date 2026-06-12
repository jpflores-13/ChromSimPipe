#!/usr/bin/env bash
# ==============================================================================
# filename:    ChromSimPipe.sh
# project:     ChromSimPipe
# description: sbatch entry point for the ChromSimPipe Snakemake pipeline.
#              Submits Snakemake as a SLURM job; Snakemake then launches all
#              downstream rules via the SLURM executor plugin.
#
# Usage:
#   sbatch ChromSimPipe.sh
#
# Requirements (run setup_data.sh first):
#   - conda envs: cohesin_sim, ctcf_extraction
#   - data/mcool/control.cool   (or let the pipeline convert it)
#   - data/ctcf_beds/*.bed       (or let the pipeline extract them)
#   - config/ChromSimConfig.yaml (edit paths before running)
#
# To run a dry-run first (shows all jobs without submitting):
#   bash ChromSimPipe.sh --dry-run
# ==============================================================================

#SBATCH --job-name=ChromSimPipe
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --partition=general

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
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
echo ""

# ── Launch pipeline ────────────────────────────────────────────────────────────
snakemake \
    --snakefile    workflows/ChromSimPipe.snakefile \
    --profile      profiles/slurm \
    --configfile   config/ChromSimConfig.yaml \
    --jobs         100 \
    --rerun-incomplete \
    --latency-wait 500 \
    --keep-going \
    ${DRY_RUN}

echo ""
echo "Pipeline complete."
