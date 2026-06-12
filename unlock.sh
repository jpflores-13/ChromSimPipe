#!/usr/bin/env bash
# ==============================================================================
# filename:    unlock.sh
# description: Release the Snakemake directory lock after a failed or
#              interrupted run. Safe to run before re-submitting SimPipe.sh.
#
# Usage:
#   bash unlock.sh
# ==============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

VENV_DIR="${REPO_ROOT}/.snakemake_venv"

if [[ -f "${VENV_DIR}/bin/activate" ]]; then
    source "${VENV_DIR}/bin/activate"
else
    module load python/3.12.4 2>/dev/null || true
fi

SNAKEFILE="workflows/ChromSimPipe.snakefile"
CONFIG="config/ChromSimConfig.yaml"
[ -f "workflows/STRSsim.snakefile" ] && SNAKEFILE="workflows/STRSsim.snakefile"
[ -f "config/SimConfig.yaml" ] && CONFIG="config/SimConfig.yaml"

snakemake \
    --snakefile  "${SNAKEFILE}" \
    --configfile "${CONFIG}" \
    --unlock

echo "Lock released."
