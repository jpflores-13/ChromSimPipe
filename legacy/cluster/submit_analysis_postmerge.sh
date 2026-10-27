#!/bin/bash
#SBATCH --job-name=analysis_postmerge
#SBATCH --output=logs/analysis_postmerge_%j.out
#SBATCH --error=logs/analysis_postmerge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=workq
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#SBATCH --mem=360G
#SBATCH --time=24:00:00

# ============================================================================
# POST-MERGE DOWNSTREAM ANALYSIS
#
# Runs the FULL post-simulation analysis pipeline on already-merged
# directories under results/polychrom_3d/merged_*/. Does NOT call
# scripts/merge_shards.py — use cluster/submit_merge_shards.sh first.
#
# This is the wrapper to use once a fresh batch of merged simulation
# directories is on disk and you want to run / re-run all downstream
# analysis without retouching the (slow) merge step:
#
#   1) sbatch cluster/submit_merge_shards.sh                # only if needed
#   2) sbatch cluster/submit_analysis_postmerge.sh          # this script
#
# What runs (matches scripts/run_analysis_all.py phases):
#   • Per-replicate: contact maps, P(s) curve + power-law fit,
#                    APA + AbLE on convergent CTCF anchors,
#                    CTCF arrow overlay, two-point MSD on probe pairs,
#                    polymer dynamics (Rg, looped fraction, dwell time).
#   • Pooled across replicates per condition: pooled contact map,
#                    pooled P(s) and derivative, pooled MSD.
#   • Cross-condition aggregates: ps_overlay_all_conditions.png,
#                    able_conserved_pairs_*, summary_table.csv,
#                    MSD α/K_α statistics.
#   • Final: contact-map figure panels via plot_contact_map_panels.py.
#
# Usage examples:
#   mkdir -p logs
#   sbatch cluster/submit_analysis_postmerge.sh
#   sbatch cluster/submit_analysis_postmerge.sh --skip-existing
#   sbatch cluster/submit_analysis_postmerge.sh --condition merged_mESC_ctcf-mESC_rep0
#   RESUME=1                sbatch cluster/submit_analysis_postmerge.sh
#   SKIP_MSD=1 SKIP_APA=1   sbatch cluster/submit_analysis_postmerge.sh
#   CLEANUP_LEGACY=1        sbatch cluster/submit_analysis_postmerge.sh
#
# Output goes to results/analysis/ as flat files (per CLAUDE.md path
# invariants); contact-map panels go to results/figures/.
# ============================================================================

set -euo pipefail

# --- SLURM-var fallbacks ------------------------------------------------------
: "${SLURM_SUBMIT_DIR:=$(pwd)}"
: "${SLURM_CPUS_PER_TASK:=32}"
: "${SLURM_MEM_PER_NODE:=}"
: "${SLURM_JOB_ID:=}"

# --- Environment --------------------------------------------------------------
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

if [ "${CONDA_DEFAULT_ENV:-}" != "polychrom" ]; then
    echo "ERROR: conda env is '${CONDA_DEFAULT_ENV:-<none>}', expected 'polychrom'." >&2
    echo "       Python in use: $(command -v python)" >&2
    exit 1
fi

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

# --- Configuration ------------------------------------------------------------
RESULTS_DIR="results/polychrom_3d"
OUTPUT_DIR="results/analysis"
HIC_DIR="data"
MCOOL_MESC="data/mcool/ES_chr3.mcool"
MCOOL_CN="data/mcool/CN_chr3.mcool"

# Conserved-CTCF BED (regenerated below unless SKIP_CONSERVED_BED=1)
CTCF_BED_MESC="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"
CTCF_BED_NEURON="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"

# Optional non-CTCF "sticky" overlay (left empty by default)
ELEMENTS_BED_MESC=""
ELEMENTS_BED_CN=""
ELEMENTS_LABEL=""

CALIBRATE_WITH="hic"
EXPT_MSD_MESC=""
EXPT_MSD_CN=""
N_BOOT=10000
STATS_SEED=0

# 16 inner workers per dir × 2 outer dirs = 32 cores fully used.
N_JOBS_PER_DIR=16

# --- Sanity check: refuse to run if no merged dirs are present ----------------
# This is what makes this wrapper "post-merge". If you genuinely want to
# analyse un-merged shards, point --results-dir somewhere else manually.
shopt -s nullglob
MERGED_DIRS=( "${RESULTS_DIR}"/merged_*/ )
shopt -u nullglob
N_MERGED=${#MERGED_DIRS[@]}
if [ "${N_MERGED}" -eq 0 ]; then
    echo "ERROR: no merged_* directories under ${RESULTS_DIR}/." >&2
    echo "       Run cluster/submit_merge_shards.sh first, or set" >&2
    echo "       RESULTS_DIR to a folder that contains merged_*/ dirs." >&2
    exit 2
fi
echo "Found ${N_MERGED} merged simulation directories under ${RESULTS_DIR}/."

# --- Print job info -----------------------------------------------------------
echo "============================================================"
echo "POST-MERGE ANALYSIS - $(date)"
echo "============================================================"
echo "Job ID:          ${SLURM_JOB_ID:-<interactive>}"
echo "Node:            $(hostname)"
echo "CPUs alloc:      ${SLURM_CPUS_PER_TASK}"
echo "Mem alloc:       ${SLURM_MEM_PER_NODE:-<n/a>} MB"
echo "Jobs per dir:    ${N_JOBS_PER_DIR}"
N_CONCURRENT=$(( SLURM_CPUS_PER_TASK / N_JOBS_PER_DIR ))
[ "${N_CONCURRENT}" -lt 1 ] && N_CONCURRENT=1
echo "Concurrent dirs: ${N_CONCURRENT}"
echo "Results dir:     ${RESULTS_DIR}"
echo "Output dir:      ${OUTPUT_DIR}"
echo "Merged dirs:     ${N_MERGED}"
echo "Resume mode:     ${RESUME:-0}  (1 = --reuse-heavy)"
echo "Extra args:      ${*:-<none>}"
echo "============================================================"

# --- Regenerate conserved CTCF BED -------------------------------------------
if [ "${SKIP_CONSERVED_BED:-0}" != "1" ] && [ -f scripts/export_conserved_ctcf_bed.py ]; then
    echo "--- Regenerating conserved CTCF BED ---"
    if python scripts/export_conserved_ctcf_bed.py; then
        echo "Conserved BED OK: ${CTCF_BED_MESC}"
    else
        echo "WARNING: export_conserved_ctcf_bed.py failed; falling back to" >&2
        echo "         per-cell-type BEDs in configs/parameters.py" >&2
        CTCF_BED_MESC=""
        CTCF_BED_NEURON=""
    fi
fi

# --- Build argument list ------------------------------------------------------
ARGS=(
    --results-dir    "${RESULTS_DIR}"
    --output-dir     "${OUTPUT_DIR}"
    --n-jobs         "${N_JOBS_PER_DIR}"
    --calibrate-with "${CALIBRATE_WITH}"
    --n-boot         "${N_BOOT}"
    --stats-seed     "${STATS_SEED}"
)

if [ -n "${HIC_DIR}" ] && [ -d "${HIC_DIR}" ] && ls "${HIC_DIR}"/*.npy &>/dev/null; then
    ARGS+=(--hic-dir "${HIC_DIR}")
fi
[ -n "${MCOOL_MESC}" ] && [ -f "${MCOOL_MESC}" ] && ARGS+=(--mcool-mesc   "${MCOOL_MESC}")
[ -n "${MCOOL_CN}"   ] && [ -f "${MCOOL_CN}"   ] && ARGS+=(--mcool-neuron "${MCOOL_CN}")
[ -n "${CTCF_BED_MESC}"   ] && [ -f "${CTCF_BED_MESC}"   ] && ARGS+=(--ctcf-bed-mesc   "${CTCF_BED_MESC}")
[ -n "${CTCF_BED_NEURON}" ] && [ -f "${CTCF_BED_NEURON}" ] && ARGS+=(--ctcf-bed-neuron "${CTCF_BED_NEURON}")
[ -n "${ELEMENTS_BED_MESC}" ] && [ -f "${ELEMENTS_BED_MESC}" ] && ARGS+=(--elements-bed-mesc   "${ELEMENTS_BED_MESC}")
[ -n "${ELEMENTS_BED_CN}"   ] && [ -f "${ELEMENTS_BED_CN}"   ] && ARGS+=(--elements-bed-neuron "${ELEMENTS_BED_CN}")
[ -n "${ELEMENTS_LABEL}" ] && ARGS+=(--elements-label "${ELEMENTS_LABEL}")
[ -n "${EXPT_MSD_MESC}" ] && [ -f "${EXPT_MSD_MESC}" ] && ARGS+=(--expt-msd-mesc   "${EXPT_MSD_MESC}")
[ -n "${EXPT_MSD_CN}"   ] && [ -f "${EXPT_MSD_CN}"   ] && ARGS+=(--expt-msd-neuron "${EXPT_MSD_CN}")

# Phase toggles (env=1 to skip)
[ "${SKIP_APA:-0}"              = "1" ] && ARGS+=(--no-apa)
[ "${SKIP_CTCF_OVERLAY:-0}"     = "1" ] && ARGS+=(--no-ctcf-overlay)
[ "${SKIP_POOL:-0}"             = "1" ] && ARGS+=(--no-pool)
[ "${SKIP_EXISTING:-0}"         = "1" ] && ARGS+=(--skip-existing)
[ "${SKIP_MSD:-0}"              = "1" ] && ARGS+=(--no-msd)
[ "${SKIP_POLYMER_DYNAMICS:-0}" = "1" ] && ARGS+=(--no-polymer-dynamics)
[ "${SKIP_MSD_STATS:-0}"        = "1" ] && ARGS+=(--no-msd-stats)
[ "${SKIP_SUMMARY:-0}"          = "1" ] && ARGS+=(--no-summary-table)
[ "${NO_PARALLEL:-0}"           = "1" ] && ARGS+=(--no-parallel)
[ "${RESUME:-0}"                = "1" ] && ARGS+=(--reuse-heavy)
[ "${CLEANUP_LEGACY:-0}"        = "1" ] && ARGS+=(--cleanup-legacy-cache)

# Forward extra CLI args verbatim
[ "$#" -gt 0 ] && ARGS+=("$@")

# --- Memory hygiene -----------------------------------------------------------
# See cluster/submit_analysis_all_32cpu.sh and CLAUDE.md "Memory landmines"
# for why these are exactly these values. Lowering MALLOC_ARENA_MAX or
# raising N_JOBS_PER_DIR < 16 has been responsible for a 360 GB cgroup OOM
# in this exact codepath; do not "tune" this without reading that section.
export MALLOC_ARENA_MAX=2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running: python scripts/run_analysis_all.py ${ARGS[*]}"
echo
python scripts/run_analysis_all.py "${ARGS[@]}"

echo
echo "============================================================"
echo "ANALYSIS DONE - $(date)"
echo "============================================================"

# --- Re-render contact-map figure panels -------------------------------------
echo
echo "--- Generating contact-map figure panels ---"
python scripts/plot_contact_map_panels.py \
    --analysis-dir "${OUTPUT_DIR}" \
    --output-dir   results/figures \
    || echo "(plot_contact_map_panels.py failed; non-fatal)"

echo
echo "ALL DONE - $(date)"
