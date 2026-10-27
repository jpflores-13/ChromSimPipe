#!/bin/bash
#SBATCH --job-name=analysis_all_16
#SBATCH --output=logs/analysis_all_16_%j.out
#SBATCH --error=logs/analysis_all_16_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# -- PARTITION -------------------------------------------------------------
# CPU-only job (no GPU). `workq` is the default CPU partition on this cluster
# (72 cores, 377 GB RAM per node, no per-user CPU QOS cap). Run
#   bash cluster/check_partitions.sh
# to confirm the exact name on your site, or swap to `standard` / `cpu`.
#
#SBATCH --partition=workq
#
# -- CPUs ------------------------------------------------------------------
# 16 cores is the right size. With SLURM_CPUS_PER_TASK=16 and
# N_JOBS_PER_DIR=16 the orchestrator derives
#   outer workers = floor(16 / 16) = 1         -> dirs run sequentially
#   inner workers = 16                         -> full per-dir parallelism
#
# --hint=nomultithread pins us to 16 *physical* cores rather than a mix of
# physical + SMT siblings. Cohesin-sim cKDTree sweeps are memory-bandwidth
# bound; SMT siblings halve per-core bandwidth for no CPU-time gain and
# inflate RSS via per-HT-thread stacks.
#
#SBATCH --cpus-per-task=16
#SBATCH --hint=nomultithread
#
# -- MEMORY ----------------------------------------------------------------
# 360 GB of 377 GB per workq node. With the 2026-04-22 Hansen-lab file-
# parallel rewrite of analysis/contact_maps.py, steady-state RSS during
# contact-map extraction is just 16 workers × ~50 MB + parent ~1 GB, so
# 360 GB is massive headroom. The previous headline of "230 GB per dir"
# reflected the pre-rewrite streaming path where the parent retained every
# pairs ndarray returned by imap_unordered and leaked ~590 KB/frame
# (job 9380631, 2026-04-22).
#
#SBATCH --mem=360G
#
# -- TIME ------------------------------------------------------------------
# 48h walltime. At ~1 h per directory (observed: 145,600 frames in ~60 min)
# a fresh 17-directory sweep finishes in ~17 h sequential. A RESUME=1 pass
# (cached contact maps + MSD, just redoing downstream figures) completes
# in under 1 h. 48 h leaves generous headroom for larger future sweeps.
#
#SBATCH --time=48:00:00

# ============================================================================
# POST-SIMULATION ANALYSIS - 16-CPU SEQUENTIAL BATCH JOB (all conditions)
#
# Runs scripts/run_analysis_all.py across EVERY merged simulation directory
# found under results/polychrom_3d/. The orchestrator auto-discovers all
# conditions (defined in configs/parameters.py) and every replicate.
#
# Output goes to results/analysis/ as a flat folder of
# {condition}_{n_blocks}blk_rep{rep}_*.{npy,png,json,...} files plus the
# pooled-replicate {condition}_*blk_pooled_*.* and the cross-condition
# aggregates. Override with --output-dir <DIR>.
#
# Faster alternative: cluster/submit_analysis_all_32cpu.sh allocates 32
# cores so two merged-dir analyses run concurrently (~9 h vs ~17 h for a
# 17-condition sweep). Use the 32-CPU variant unless queue pressure is
# tight. Both write to the same flat results/analysis/ folder.
#
# Pre-step: the conserved mESC-intersect-CN oriented CTCF BED is
# auto-regenerated at the top of the job via scripts/export_conserved_ctcf_bed.py
# so downstream outputs (APA pileup, contact-map arrows, MSD probe pair)
# always reflect the latest per-cell-type BED edits. Skip with
#   SKIP_CONSERVED_BED=1 sbatch cluster/submit_analysis_all_16cpu.sh
#
# For each merged directory the five core analysis modules are invoked:
#     analysis/contact_maps.py
#     analysis/ps_curve.py
#     analysis/absolute_quant.py
#     analysis/ctcf_plotting.py
#     analysis/experimental_compare.py
# plus the auxiliary modules triggered by the orchestrator:
#     analysis/msd_two_point.py, polymer_dynamics.py, calibration.py, lef_lifetimes.py
#
# After all per-replicate work finishes the orchestrator performs:
#   - pooled-replicate analysis (Hansen-lab style concat across replicates),
#   - Phase-3 between-condition MSD statistics (forest plots + p-values),
#   - a combined summary table across all conditions.
#
# Usage:
#   mkdir -p logs
#   sbatch cluster/submit_analysis_all_16cpu.sh
#   sbatch cluster/submit_analysis_all_16cpu.sh --skip-existing
#   sbatch cluster/submit_analysis_all_16cpu.sh --condition <merged_dir_name>
#   SKIP_MSD=1 SKIP_POLYMER_DYNAMICS=1 sbatch cluster/submit_analysis_all_16cpu.sh
#   RESUME=1 sbatch cluster/submit_analysis_all_16cpu.sh    # 2nd pass: reuse
#                                                           # contact maps + MSD,
#                                                           # redo APA / CTCF /
#                                                           # summaries only
#   SKIP_CONSERVED_BED=1 sbatch cluster/submit_analysis_all_16cpu.sh  # skip pre-step
#   CLEANUP_LEGACY=1 sbatch cluster/submit_analysis_all_16cpu.sh      # also delete
#                                                          # leftover
#                                                          # <sim_dir>/analysis/
#                                                          # directories from
#                                                          # the old layout
#
# Every flag documented in the "CLI flags for run_analysis_all.py" table in
# README.md is supported as a passthrough.
# ============================================================================

set -euo pipefail

# --- SLURM-var fallbacks ------------------------------------------------------
: "${SLURM_SUBMIT_DIR:=$(pwd)}"
: "${SLURM_CPUS_PER_TASK:=16}"
: "${SLURM_MEM_PER_NODE:=}"
: "${SLURM_JOB_ID:=}"

# --- Environment --------------------------------------------------------------
# IMPORTANT: in a non-interactive sbatch shell `conda activate` only works
# AFTER `eval "$(conda shell.bash hook)"`. Without the hook, the `activate`
# call silently falls through to base - the job then runs on base's Python
# (3.10) instead of polychrom's (3.12), producing confusing
#   ImportError: cannot import name '...' from '...' (unknown location)
# errors deep in the analysis modules.
export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate polychrom

# Fail loud if activation did not land us in polychrom.
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
MCOOL_MESC="data/mcool/ES_chr3.mcool"  # Option A chr3-only default
MCOOL_CN="data/mcool/CN_chr3.mcool"

# --- Oriented CTCF BED overrides ---------------------------------------------
# By default the orchestrator uses the per-cell-type BEDs declared in
# configs/parameters.py (CTCF_BED_MESC / CTCF_BED_NEURON). Setting BOTH of
# the variables below to the SAME "conserved" BED makes every downstream
# analysis (APA pileup, contact-map CTCF arrows, and - indirectly via the
# intersection rule already in parameters.py - the MSD probe pair) focus on
# the mESC intersect CN conserved CTCF sites.
#
# The conserved BED file pointed to below is (re)generated automatically at
# the top of this job via scripts/export_conserved_ctcf_bed.py.
# Leave either variable empty ("") to fall back to that cell type's default
# BED (the auto-discovered configs.parameters value).
CTCF_BED_MESC="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"
CTCF_BED_NEURON="data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"

# --- Optional non-CTCF "sticky" overlay --------------------------------------
ELEMENTS_BED_MESC=""
ELEMENTS_BED_CN=""
ELEMENTS_LABEL=""

CALIBRATE_WITH="hic"                   # hic | msd | none
EXPT_MSD_MESC=""
EXPT_MSD_CN=""
N_BOOT=10000
STATS_SEED=0

# --- Inner-pool parallelism ---------------------------------------------------
# Number of workers in the inner multiprocessing.Pool used by
# analysis/contact_maps.py for cKDTree contact detection (and by a handful
# of other per-dir modules that honour the same knob).
#
# Keep equal to SLURM_CPUS_PER_TASK so the orchestrator's outer worker
# count collapses to 1:
#   n_dir_workers = floor(16 / 16) = 1
# i.e. directories run STRICTLY sequentially, with 16 inner Pool workers
# each. Same 16 cores used, different layout.
#
# Why this matters (2026-04-22 post-mortem of 9360298/9374903/9375923/9377438):
# Phase 1 parallelises directories via a *ThreadPoolExecutor*
# (scripts/run_analysis_all.py:1575). Threads share one Python heap, so
# with N_JOBS_PER_DIR=4 four threads concurrently allocate numpy
# temporaries in the same process while each also holds a
# multiprocessing.Pool(4, maxtasksperchild=500). Every worker recycle is
# a fork() that inherits that heap via CoW; over ~95k tasks per dir,
# glibc arena fragmentation + un-shared CoW pages drove RSS past
# --mem=360G and SLURM cgroup-OOM-killed the step. Reverting to a single
# outer thread eliminates the concurrent-allocator thrash on the parent
# heap; each new fork inherits a much smaller, quieter RSS.
#
# Sequential dirs are fine walltime-wise: ~1 h/dir × 17 dirs ≈ 17 h,
# well inside the 48 h budget.
N_JOBS_PER_DIR=16

# --- Feature toggles (set to 1 via environment to disable) --------------------
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
echo "POST-SIM ANALYSIS (16-CPU sequential, all conditions) - $(date)"
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
# CTCF_BED_MESC and CTCF_BED_NEURON above point to the mESC intersect CN
# "conserved" oriented BED. That file is derived from the two per-cell-type
# BEDs declared in configs/parameters.py. We rebuild it every time the
# analysis is submitted, so any upstream edit to either source BED is
# automatically reflected in every downstream output (APA pileup,
# contact-map CTCF arrows, MSD probe-pair selection).
#
# Skip with   SKIP_CONSERVED_BED=1 sbatch ...   to keep the existing file
# (useful during rapid iteration on downstream plotting only).
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

# Append any positional args the user passed to sbatch (e.g. --skip-existing,
# --condition <name>). Quoted "$@" preserves args that contain spaces —
# the previous `ARGS+=( ${EXTRA_ARGS} )` split on IFS and broke such args.
if [ "$#" -gt 0 ]; then
    ARGS+=("$@")
fi

# --- Run ----------------------------------------------------------------------
# Cap glibc malloc arenas. Default = 8 × ncores = 128 on workq, each arena
# can waste up to ~64 MB in per-thread fragmentation caches → ~8 GB of
# waste per process. cKDTree builds lots of small allocations inside each
# worker, which is exactly the pattern ptmalloc fragments on. Two arenas
# is the recommended value for numpy/scipy parallel workloads and keeps
# the 16 inner workers' RSS growth bounded over long (2.8M-task) runs.
export MALLOC_ARENA_MAX=2

# Pin BLAS / OpenMP thread pools to 1 thread per worker process. Without
# this, each of the 16 multiprocessing.Pool workers would see
# CPU_COUNT=16 and spawn 16 internal numpy/scipy BLAS threads → 256 OS
# threads fighting over 16 hardware cores. That inflates RSS via
# per-thread stacks and trashes memory bandwidth for cKDTree sweeps.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Running: MALLOC_ARENA_MAX=${MALLOC_ARENA_MAX} OMP_NUM_THREADS=${OMP_NUM_THREADS} python scripts/run_analysis_all.py ${ARGS[*]}"
echo ""
python scripts/run_analysis_all.py "${ARGS[@]}"

# Regenerate the contact-map figure panels (per-rep grids, pooled
# merged_reps, all-conditions overview) into results/figures/.
echo ""
echo "--- Generating contact-map figure panels ---"
python scripts/plot_contact_map_panels.py \
    --analysis-dir "${OUTPUT_DIR}" \
    --output-dir results/figures || echo "(plot_contact_map_panels.py failed; non-fatal)"

echo ""
echo "============================================================"
echo "DONE - $(date)"
echo "============================================================"
