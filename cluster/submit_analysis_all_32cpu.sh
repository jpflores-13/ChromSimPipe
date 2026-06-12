#!/bin/bash
#SBATCH --job-name=analysis_all_32
#SBATCH --output=logs/analysis_all_32_%j.out
#SBATCH --error=logs/analysis_all_32_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# -- PARTITION -------------------------------------------------------------
# CPU-only job. workq = 72 cores, 377 GB RAM per node, no per-user CPU QOS
# cap. Run `bash cluster/check_partitions.sh` if unsure of the partition
# name on this site.
#
#SBATCH --partition=workq
#
# -- CPUs ------------------------------------------------------------------
# 32 cores. With SLURM_CPUS_PER_TASK=32 and N_JOBS_PER_DIR=16 the
# orchestrator derives:
#   outer workers = floor(32 / 16) = 2          -> TWO dirs run in parallel
#   inner workers = 16 each                     -> full per-dir parallelism
#
# This roughly halves wallclock vs. the 16-CPU sequential variant
# (cluster/submit_analysis_all_16cpu.sh) for a 17-condition sweep.
# Net cost: 32 cores allocated for ~9 h instead of 16 cores for ~17 h —
# same CPU-h budget, ~2x faster turnaround.
#
# Why two parallel dirs is safe NOW (was not before 2026-04-22):
# The pre-2026-04-22 contact-map extractor accumulated every frame's
# pairs ndarray on the parent heap (~590 KB/frame leak via glibc
# fragmentation). With ~95k frames per dir, two concurrent dirs would
# pile up >100 GB of un-shared CoW pages and trip --mem cgroup OOM.
# After the Hansen-lab file-parallel rewrite of analysis/contact_maps.py
# (commit 9380631), steady-state RSS during contact-map extraction is
# only 16 workers × ~50 MB + parent ~1 GB = ~1.8 GB per dir. Two dirs
# = ~3.6 GB. 360 GB --mem is massive headroom.
#
# --hint=nomultithread pins us to 32 *physical* cores. Cohesin-sim
# cKDTree sweeps are memory-bandwidth bound; SMT siblings halve
# per-core bandwidth for no CPU-time gain.
#
#SBATCH --cpus-per-task=32
#SBATCH --hint=nomultithread
#
# -- MEMORY ----------------------------------------------------------------
# 360 GB of 377 GB per workq node. Same allocation as the 16-CPU variant;
# steady-state usage with 2 parallel dirs is ~4 GB total.
#
#SBATCH --mem=360G
#
# -- TIME ------------------------------------------------------------------
# 24h walltime. With 2 parallel dirs at ~1 h each, a 17-dir sweep
# completes in ~9 h. RESUME=1 (cached contact maps + MSD) finishes in
# under 1 h. 24 h leaves comfortable headroom.
#
#SBATCH --time=24:00:00

# ============================================================================
# POST-SIMULATION ANALYSIS - 32-CPU 2x-PARALLEL BATCH JOB (all conditions)
#
# Same script as submit_analysis_all_16cpu.sh but allocates 32 cores so the
# orchestrator runs two merged-directory analyses in parallel. Use this
# when you want to halve wallclock at no extra CPU-h cost. Use the 16-CPU
# variant when the cluster is full and queue time matters more than wall
# time.
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_analysis_all_32cpu.sh
#   sbatch cluster/submit_analysis_all_32cpu.sh --skip-existing
#   sbatch cluster/submit_analysis_all_32cpu.sh --condition <merged_dir>
#   RESUME=1 sbatch cluster/submit_analysis_all_32cpu.sh
#   CLEANUP_LEGACY=1 sbatch cluster/submit_analysis_all_32cpu.sh
#
# Output goes to results/analysis/ as a flat folder of
# {condition}_{n_blocks}blk_rep{rep}_*.{npy,png,json,...} files plus the
# pooled-replicate {condition}_*blk_pooled_*.* and the cross-condition
# aggregates (ps_overlay_all_conditions.png, summary_table.csv, ...).
# Override the destination with --output-dir <DIR>.
# ============================================================================

set -euo pipefail

# --- SLURM-var fallbacks ------------------------------------------------------
: "${SLURM_SUBMIT_DIR:=$(pwd)}"
: "${SLURM_CPUS_PER_TASK:=32}"
: "${SLURM_MEM_PER_NODE:=}"
: "${SLURM_JOB_ID:=}"

# --- Environment --------------------------------------------------------------
# IMPORTANT: in a non-interactive sbatch shell `conda activate` only works
# AFTER `eval "$(conda shell.bash hook)"`. Without the hook, the `activate`
# call silently falls through to base.
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

if [ "${CONDA_DEFAULT_ENV:-}" != "polychrom" ]; then
    echo "ERROR: conda env is '${CONDA_DEFAULT_ENV:-<none>}', expected 'polychrom'." >&2
    echo "       Python in use: $(command -v python)" >&2
    exit 1
fi
echo "Python:        $(command -v python)"
echo "Python ver:    $(python --version 2>&1)"
echo "Conda env:     ${CONDA_DEFAULT_ENV}"

# --- Configuration ------------------------------------------------------------
RESULTS_DIR="results/polychrom_3d"
OUTPUT_DIR="results/analysis"          # flat output folder (also the default)
HIC_DIR="data"                         # legacy hic_*_Sox2.npy; "" disables
MCOOL_MESC="data/mcool/ES_chr3.mcool"
MCOOL_CN="data/mcool/CN_chr3.mcool"

# --- Oriented CTCF BED overrides ---------------------------------------------
CTCF_BED_MESC="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"
CTCF_BED_NEURON="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"

# --- Optional non-CTCF "sticky" overlay --------------------------------------
ELEMENTS_BED_MESC=""
ELEMENTS_BED_CN=""
ELEMENTS_LABEL=""

CALIBRATE_WITH="hic"
EXPT_MSD_MESC=""
EXPT_MSD_CN=""
N_BOOT=10000
STATS_SEED=0

# --- Inner-pool parallelism ---------------------------------------------------
# 16 inner workers per dir × 2 outer dirs = 32 cores fully used.
N_JOBS_PER_DIR=16

# --- Feature toggles (env=1 to enable) ---------------------------------------
# SKIP_APA=1 SKIP_CTCF_OVERLAY=1 SKIP_POOL=1 SKIP_EXISTING=1
# SKIP_MSD=1 SKIP_POLYMER_DYNAMICS=1 SKIP_MSD_STATS=1 NO_PARALLEL=1
# SKIP_CONSERVED_BED=1   -> do not regenerate the conserved CTCF BED
# CLEANUP_LEGACY=1       -> delete <sim_dir>/analysis/ after each rep
#
# RESUME=1  ->  pass --reuse-heavy to run_analysis_all.py.

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

# --- Print job info -----------------------------------------------------------
echo "============================================================"
echo "POST-SIM ANALYSIS (32-CPU 2x-parallel, all conditions) - $(date)"
echo "============================================================"
echo "Job ID:          ${SLURM_JOB_ID:-<interactive>}"
echo "Node:            $(hostname)"
echo "CPUs alloc:      ${SLURM_CPUS_PER_TASK}"
echo "Mem alloc:       ${SLURM_MEM_PER_NODE:-<n/a>} MB"
echo "Jobs per dir:    ${N_JOBS_PER_DIR}"
N_CONCURRENT=$(( SLURM_CPUS_PER_TASK / N_JOBS_PER_DIR ))
[ "${N_CONCURRENT}" -lt 1 ] && N_CONCURRENT=1
echo "Concurrent dirs: ${N_CONCURRENT}  (= ${SLURM_CPUS_PER_TASK} / ${N_JOBS_PER_DIR})"
echo "Results dir:     ${RESULTS_DIR}"
echo "Output dir:      ${OUTPUT_DIR}"
echo "Hi-C dir:        ${HIC_DIR:-<none>}"
echo "mcool mESC:      ${MCOOL_MESC:-<none>}"
echo "mcool CN:        ${MCOOL_CN:-<none>}"
echo "Calibrate with:  ${CALIBRATE_WITH}"
echo "Resume mode:     ${RESUME:-0}  (1 = reuse heavy cached outputs)"
echo "Extra args:      ${*:-<none>}"
echo "============================================================"
echo ""

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
    echo ""
fi

# --- Build argument list ------------------------------------------------------
ARGS=(
    --results-dir "${RESULTS_DIR}"
    --output-dir  "${OUTPUT_DIR}"
    --n-jobs      "${N_JOBS_PER_DIR}"
    --calibrate-with "${CALIBRATE_WITH}"
    --n-boot      "${N_BOOT}"
    --stats-seed  "${STATS_SEED}"
)

if [ -n "${HIC_DIR}" ] && [ -d "${HIC_DIR}" ]; then
    if ls "${HIC_DIR}"/*.npy &>/dev/null; then
        ARGS+=(--hic-dir "${HIC_DIR}")
    else
        echo "NOTE: no .npy files in ${HIC_DIR} - skipping legacy .npy comparison"
    fi
fi

if [ -n "${MCOOL_MESC}" ] && [ -f "${MCOOL_MESC}" ]; then
    ARGS+=(--mcool-mesc "${MCOOL_MESC}")
fi
if [ -n "${MCOOL_CN}" ] && [ -f "${MCOOL_CN}" ]; then
    ARGS+=(--mcool-neuron "${MCOOL_CN}")
fi

if [ -n "${CTCF_BED_MESC}" ] && [ -f "${CTCF_BED_MESC}" ]; then
    ARGS+=(--ctcf-bed-mesc "${CTCF_BED_MESC}")
fi
if [ -n "${CTCF_BED_NEURON}" ] && [ -f "${CTCF_BED_NEURON}" ]; then
    ARGS+=(--ctcf-bed-neuron "${CTCF_BED_NEURON}")
fi

if [ -n "${ELEMENTS_BED_MESC}" ] && [ -f "${ELEMENTS_BED_MESC}" ]; then
    ARGS+=(--elements-bed-mesc "${ELEMENTS_BED_MESC}")
fi
if [ -n "${ELEMENTS_BED_CN}" ] && [ -f "${ELEMENTS_BED_CN}" ]; then
    ARGS+=(--elements-bed-neuron "${ELEMENTS_BED_CN}")
fi
if [ -n "${ELEMENTS_LABEL}" ]; then
    ARGS+=(--elements-label "${ELEMENTS_LABEL}")
fi

if [ -n "${EXPT_MSD_MESC}" ] && [ -f "${EXPT_MSD_MESC}" ]; then
    ARGS+=(--expt-msd-mesc "${EXPT_MSD_MESC}")
fi
if [ -n "${EXPT_MSD_CN}" ] && [ -f "${EXPT_MSD_CN}" ]; then
    ARGS+=(--expt-msd-neuron "${EXPT_MSD_CN}")
fi

[ "${SKIP_APA:-0}"              = "1" ] && ARGS+=(--no-apa)
[ "${SKIP_CTCF_OVERLAY:-0}"     = "1" ] && ARGS+=(--no-ctcf-overlay)
[ "${SKIP_POOL:-0}"             = "1" ] && ARGS+=(--no-pool)
[ "${SKIP_EXISTING:-0}"         = "1" ] && ARGS+=(--skip-existing)
[ "${SKIP_MSD:-0}"              = "1" ] && ARGS+=(--no-msd)
[ "${SKIP_POLYMER_DYNAMICS:-0}" = "1" ] && ARGS+=(--no-polymer-dynamics)
[ "${SKIP_MSD_STATS:-0}"        = "1" ] && ARGS+=(--no-msd-stats)
[ "${NO_PARALLEL:-0}"           = "1" ] && ARGS+=(--no-parallel)
[ "${RESUME:-0}"                = "1" ] && ARGS+=(--reuse-heavy)
[ "${CLEANUP_LEGACY:-0}"        = "1" ] && ARGS+=(--cleanup-legacy-cache)

if [ "$#" -gt 0 ]; then
    ARGS+=("$@")
fi

# --- Run ----------------------------------------------------------------------
# Cap glibc malloc arenas. Default = 8 × ncores = 256 on 32-core. Two
# arenas is the recommended value for numpy/scipy parallel workloads.
export MALLOC_ARENA_MAX=2

# Pin BLAS / OpenMP thread pools to 1 thread per worker process.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running: MALLOC_ARENA_MAX=${MALLOC_ARENA_MAX} OMP_NUM_THREADS=${OMP_NUM_THREADS} python scripts/run_analysis_all.py ${ARGS[*]}"
echo ""
python scripts/run_analysis_all.py "${ARGS[@]}"

echo ""
echo "============================================================"
echo "DONE - $(date)"
echo "============================================================"

# Optional follow-up: regenerate the contact-map figure panels.
echo ""
echo "--- Generating contact-map figure panels ---"
python scripts/plot_contact_map_panels.py \
    --analysis-dir "${OUTPUT_DIR}" \
    --output-dir results/figures || echo "(plot_contact_map_panels.py failed; non-fatal)"
