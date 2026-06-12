#!/bin/bash
#SBATCH --job-name=contact_maps
#SBATCH --output=logs/contact_maps_%j.out
#SBATCH --error=logs/contact_maps_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# ── PARTITION ──────────────────────────────────────────────────────────────
# No GPU needed. Run check_partitions.sh first to find the right partition.
# Common names: cpu, batch, standard, short, long, compute, highmem
#
# workq is the default CPU partition: 72 cores, 377 GB RAM, no CPU QOS cap.
# Other CPU options: demux (72 cores), interactive (72 cores), neuroimaging (64 cores).
# Avoid cuda — that's the GPU partition.
#SBATCH --partition=workq
#
# ── CPUs ───────────────────────────────────────────────────────────────────
# workq nodes have 72 cores; no per-user CPU cap on this account.
# With N_JOBS_PER_DIR=4 (below): floor(72/4) = 18 directories run concurrently.
#
#SBATCH --cpus-per-task=72
#
# ── MEMORY ─────────────────────────────────────────────────────────────────
# 18 concurrent dirs × 7 GB/dir + 4 GB headroom = 130 GB.
# workq nodes have 377 GB available, so this is well within limits.
#
#SBATCH --mem=130G
#
# ── TIME ───────────────────────────────────────────────────────────────────
# 48 directories × ~40 000 frames each: allow 12-24h to be safe.
# Reduce to 4:00:00 if using --skip-existing on a mostly-done run.
#
#SBATCH --time=24:00:00

# ============================================================================
# CONTACT MAP ANALYSIS — CPU BATCH JOB
#
# Runs run_analysis_all.py across all merged simulation directories.
# No GPU required. Parallelises at the directory level (one Pool of workers
# across directories) and within each directory (cKDTree contact detection).
#
# BEFORE SUBMITTING:
#   1. Run bash cluster/check_partitions.sh  to find the right partition and
#      maximum CPUs/memory available to you.
#   2. Edit --partition, --cpus-per-task, and --mem above accordingly.
#   3. Adjust N_JOBS_PER_DIR below if needed (default 4 is conservative;
#      raise to 8 on nodes with many cores to increase per-dir throughput).
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_analysis_cpu.sh
#
# To only re-run failed / new directories:
#   sbatch cluster/submit_analysis_cpu.sh --skip-existing   # passes arg through
# ============================================================================

set -euo pipefail

# --- Environment ---
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# --- Configuration ---
RESULTS_DIR="results/polychrom_3d"
HIC_DIR="data"                  # directory with legacy hic_*_Sox2.npy files.
                                # set to "" to skip the .npy-based comparison.

# --- Experimental Hi-C via mcool (preferred over .npy when available) ---
# Leave empty to skip. Paths are resolved from the repo root.
#
# Default filenames match the Bonev chr3-only pipeline in
# scripts/download_bonev_hic.sh (fast Option A), which emits
#   data/mcool/ES_chr3.mcool  and  data/mcool/CN_chr3.mcool
# If you ran Option B (full genome-wide FASTQ → mcool pipeline instead),
# the outputs are data/mcool/ES.mcool / CN.mcool — edit these two lines.
MCOOL_MESC="data/mcool/ES_chr3.mcool"
MCOOL_CN="data/mcool/CN_chr3.mcool"

# --- Optional "sticky elements" tracks (enhancers/promoters/etc) ---
# Leave empty to omit the overlay. These BEDs are drawn as blue squares
# on the 1D tracks above and to the left of the contact map.
ELEMENTS_BED_MESC=""
ELEMENTS_BED_CN=""
ELEMENTS_LABEL="enhancers/promoters"

# --- Feature toggles (uncomment to disable) ---
# SKIP_APA=1            # skip absolute-loop-quant APA pileups
# SKIP_CTCF_OVERLAY=1   # skip the CTCF-aligned BED + figure
# SKIP_POOL=1           # skip the pooled-replicate step
# SKIP_EXISTING=1       # re-use directories that already have sim_contact_map.npy

# How many cores to use for contact detection WITHIN each directory.
# The script automatically computes how many directories run concurrently:
#   n_concurrent = floor(SLURM_CPUS_PER_TASK / N_JOBS_PER_DIR)
N_JOBS_PER_DIR=4

# Pass any extra args from the command line through (e.g. --skip-existing)
EXTRA_ARGS="${@}"

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

# --- Print job info ---
echo "============================================================"
echo "CONTACT MAP ANALYSIS — $(date)"
echo "============================================================"
echo "Node:          $(hostname)"
echo "CPUs alloc:    ${SLURM_CPUS_PER_TASK}"
echo "Jobs per dir:  ${N_JOBS_PER_DIR}"
N_CONCURRENT=$(( SLURM_CPUS_PER_TASK / N_JOBS_PER_DIR ))
echo "Concurrent dirs: ${N_CONCURRENT}  (= ${SLURM_CPUS_PER_TASK} / ${N_JOBS_PER_DIR})"
echo "Results dir:   ${RESULTS_DIR}"
echo "Hi-C dir:      ${HIC_DIR:-<none, skipping comparison>}"
echo "Extra args:    ${EXTRA_ARGS:-<none>}"
echo "============================================================"
echo ""

# --- Build argument list ---
ARGS=(
    --results-dir "${RESULTS_DIR}"
    --n-jobs      "${N_JOBS_PER_DIR}"
)

if [ -n "${HIC_DIR}" ] && [ -d "${HIC_DIR}" ]; then
    # Only add --hic-dir if at least one .npy matrix exists there
    if ls "${HIC_DIR}"/*.npy &>/dev/null; then
        ARGS+=(--hic-dir "${HIC_DIR}")
    else
        echo "NOTE: no .npy files found in ${HIC_DIR} — skipping .npy Hi-C comparison"
    fi
fi

# mcool paths → --mcool-mesc / --mcool-neuron (takes precedence over .npy)
if [ -n "${MCOOL_MESC}" ] && [ -f "${MCOOL_MESC}" ]; then
    ARGS+=(--mcool-mesc "${MCOOL_MESC}")
fi
if [ -n "${MCOOL_CN}" ] && [ -f "${MCOOL_CN}" ]; then
    ARGS+=(--mcool-neuron "${MCOOL_CN}")
fi

# Sticky-element BEDs → --elements-bed-mesc / --elements-bed-neuron
if [ -n "${ELEMENTS_BED_MESC}" ] && [ -f "${ELEMENTS_BED_MESC}" ]; then
    ARGS+=(--elements-bed-mesc "${ELEMENTS_BED_MESC}")
fi
if [ -n "${ELEMENTS_BED_CN}" ] && [ -f "${ELEMENTS_BED_CN}" ]; then
    ARGS+=(--elements-bed-neuron "${ELEMENTS_BED_CN}")
fi
if [ -n "${ELEMENTS_LABEL}" ]; then
    ARGS+=(--elements-label "${ELEMENTS_LABEL}")
fi

# Feature toggles
[ "${SKIP_APA:-0}"          = "1" ] && ARGS+=(--no-apa)
[ "${SKIP_CTCF_OVERLAY:-0}" = "1" ] && ARGS+=(--no-ctcf-overlay)
[ "${SKIP_POOL:-0}"         = "1" ] && ARGS+=(--no-pool)
[ "${SKIP_EXISTING:-0}"     = "1" ] && ARGS+=(--skip-existing)

# Append any command-line passthrough args (e.g. --skip-existing)
if [ -n "${EXTRA_ARGS}" ]; then
    ARGS+=( ${EXTRA_ARGS} )
fi

# --- Run ---
echo "Running: python scripts/run_analysis_all.py ${ARGS[*]}"
echo ""
python scripts/run_analysis_all.py "${ARGS[@]}"

echo ""
echo "============================================================"
echo "DONE — $(date)"
echo "============================================================"
