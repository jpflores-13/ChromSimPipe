#!/bin/bash
#SBATCH --job-name=catchup
#SBATCH --output=logs/catchup_%A_%a.out
#SBATCH --error=logs/catchup_%A_%a.err
#SBATCH --partition=cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --signal=B:USR1@90
#SBATCH --mail-user=gabriele.michele@hsr.it
#SBATCH --mail-type=END,FAIL,TIME_LIMIT

# ============================================================================
# submit_catch_up.sh — idempotent multi-round top-up of polychrom 3D reps
#
# Constraints encoded
#   * 2 concurrent GPUs max on the cluster  → SLURM array throttled with %2
#   * 48 h walltime per task                → matches HSR cuda partition cap
#   * Each task ALWAYS passes --resume --n-blocks ${BLOCKS_PER_REP:-70000}
#     to scripts/run_simulation.py. UNIT NOTE: --n-blocks counts polychrom
#     CONFORMATIONS (== sim.do_block calls == saved frames), NOT .h5 files;
#     polychrom HDF5Reporter packs 100 conformations per blocks_*.h5
#     (max_data_length=100), so 70 000 conformations → 700 .h5 files on
#     disk. Conflating these is how the 2026-04-30 catchup landed with
#     7 .h5 files per rep instead of 700; the var is now 70 000 by default.
#   * Resume via HDF5Reporter.continue_trajectory() picks up from the last
#     saved block. Combined with the SIGUSR1 walltime trap below, a rep
#     that hits walltime self-schedules a follow-up task with the same
#     (condition, replicate) and the new task continues the same blocks_*.h5
#     trajectory — so a rep is never abandoned just because 48 h wasn't enough.
#   * --signal=B:USR1@90 makes SLURM send SIGUSR1 to the batch script 90s
#     before the wall, giving the trap time to sbatch the resume task.
#   * --mail-user / --mail-type send email on END / FAIL / TIME_LIMIT for
#     each array task — mute by overriding to "" if the volume is too high.
#
# The pool aggregator (analysis/run_analysis_all.py:discover_merged_dirs) finds
# every non-shard rep dir with blocks_*.h5, so each fresh rep number
# accumulates into the pool just like the resumed ones.
#
# Usage
# -----
#   # Plan + submit one round (called from project root):
#   bash cluster/submit_catch_up.sh
#
#   # Dry-run — print the plan, don't sbatch:
#   PLAN_ONLY=1 bash cluster/submit_catch_up.sh
#
#   # Defaults: TARGET_FRAMES_PER_CONDITION=200000, MIN_REPS_PER_CONDITION=3,
#   # MISSING_ONLY=0, INCLUDE_NEW_CONDITIONS=1. So a no-arg invocation tops
#   # every known condition up to >=3 reps AND >=200 000 pooled frames.
#   bash cluster/submit_catch_up.sh
#
#   # Override target frames per condition (e.g. shrink the depth budget):
#   TARGET_FRAMES_PER_CONDITION=150000 bash cluster/submit_catch_up.sh
#
#   # Skip the 3 conditions that have never been run (CN_high_density,
#   # CN_long_res_high_dens, CN_long_res_low_dens):
#   INCLUDE_NEW_CONDITIONS=0 bash cluster/submit_catch_up.sh
#
#   # ONLY top up the conditions that currently have NO data:
#   MISSING_ONLY=1 bash cluster/submit_catch_up.sh
#
#   # Override the rep floor (e.g. push to 5 if you want tighter MSD CIs):
#   MIN_REPS_PER_CONDITION=5 bash cluster/submit_catch_up.sh
#
# Iteration
# ---------
#   Each invocation = one round:
#     1. PLANNER (no SLURM_ARRAY_TASK_ID): scans results/polychrom_3d/, counts
#        current frames per condition, schedules new reps to reach TARGET,
#        sbatch's WORKER as a job array %2.
#     2. WORKER (with SLURM_ARRAY_TASK_ID): reads its (condition, replicate)
#        pair from the queue file and runs scripts/run_simulation.py on 1 GPU.
#
#   After the array finishes (or hits walltime), invoke again — the planner
#   now sees the freshly-written reps and only schedules what's still missing.
#
# After everything reaches the target, refresh analysis with:
#   rm -f results/analysis/*_pooled_*
#   rm -f results/analysis/*_rep*_contact_map.npy
#   RESUME=1 sbatch cluster/submit_analysis_all_32cpu.sh
# ============================================================================

set -euo pipefail

# ---------- Tunables ----------
# Defaults reflect the agreed policy (see CLAUDE.md "Statistical-depth policy"):
# every condition pooled to ≥ 200 000 frames AND with ≥ 3 reps so that per-rep
# MSD alpha CIs are comparable across conditions. Each rep gets a deterministic
# seed via run_simulation.py:rep_seed = 42 + replicate * 1000 (covers numpy,
# the LEFSimulator, and the OpenMM Langevin integrator).
TARGET_FRAMES_PER_CONDITION="${TARGET_FRAMES_PER_CONDITION:-200000}"
EXPECTED_BLOCKS_PER_REP="${EXPECTED_BLOCKS_PER_REP:-1500}"
# Block budget per rep, passed through to run_simulation.py:--n-blocks.
# UNIT: this is polychrom blocks = conformations = sim.do_block iterations.
# polychrom HDF5Reporter packs 100 conformations per blocks_*.h5 file
# (max_data_length=100), so 70 000 here = 700 .h5 files on disk =
# 70 000 saved frames per rep. ~12 h on a V100 (well within 48 h walltime).
# 3 reps × 70 000 frames = 210 000 frames per condition, just above the
# 200 000 target because comparable_merge symlinks whole .h5 files
# (smallest unit it can pick).
#
# DO NOT confuse with EXPECTED_BLOCKS_PER_REP above: that one is in .h5
# files (which is multiplied by FRAMES_PER_BLOCK=100 when computing the
# planner's frame-deficit math). They live in different units because
# run_simulation.py's --n-blocks API and analysis's "blk" filename
# convention both count conformations, while planner deficits are easier
# to reason about in files.
BLOCKS_PER_REP="${BLOCKS_PER_REP:-70000}"
# Auto-chain comparable_merge + analysis after the array drains.
# Set CHAIN_FOLLOW_UPS=0 to disable (e.g. for ad-hoc partial runs).
CHAIN_FOLLOW_UPS="${CHAIN_FOLLOW_UPS:-1}"
INCLUDE_NEW_CONDITIONS="${INCLUDE_NEW_CONDITIONS:-1}"
MIN_REPS_PER_CONDITION="${MIN_REPS_PER_CONDITION:-3}"   # floor on rep count per condition (matters for MSD CIs)
MISSING_ONLY="${MISSING_ONLY:-0}"                       # set 1 to skip conditions that already have data
RESULTS_BASE="results/polychrom_3d"
QUEUE_FILE="cluster/_catch_up_queue.txt"
FRAMES_PER_BLOCK=100   # polychrom HDF5Reporter max_data_length

# Map: on-disk dir stem → CLI condition name. The disk stem is built by
# scripts/run_simulation.py:163 as f"{params['name']}_ctcf-{ctcf_type}", so
# the CLI name is the SIMULATION_CONDITIONS entry that yields that combo.
declare -A DISK_TO_CONDITION=(
    ["mESC_ctcf-mESC"]="mESC_ctrl"
    ["mESC_ctcf-neuron"]="mESC_params_neuron_ctcf"
    ["CN_baseline_ctcf-neuron"]="CN_baseline_neuron_ctcf"
    ["CN_long_residency_ctcf-neuron"]="CN_long_residency_neuron_ctcf"
    ["CN_very_long_residency_ctcf-neuron"]="CN_very_long_residency_neuron_ctcf"
    ["CN_high_density_ctcf-neuron"]="CN_high_density_neuron_ctcf"
    ["CN_long_res_high_dens_ctcf-neuron"]="CN_long_res_high_dens_neuron_ctcf"
    ["CN_long_res_low_dens_ctcf-neuron"]="CN_long_res_low_dens_neuron_ctcf"
)

# ---------- Worker mode ----------
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    cd "${SLURM_SUBMIT_DIR}"
    export PATH=/opt/common/tools/ric.tiget/mambaforge/bin/:$PATH
    eval "$(conda shell.bash hook)"
    conda activate polychrom

    # A task can be invoked two ways:
    # 1. As an element of the planner-submitted array → reads its line from
    #    QUEUE_FILE by SLURM_ARRAY_TASK_ID.
    # 2. As a self-resubmitted resume task → the resume trap below sbatch's
    #    a 1-element array and exports RESUME_COND and RESUME_REP, which
    #    take precedence over the queue lookup.
    if [ -n "${RESUME_COND:-}" ] && [ -n "${RESUME_REP:-}" ]; then
        cond_cli="${RESUME_COND}"
        rep="${RESUME_REP}"
        echo "[resume task] using exported (cond, rep) = (${cond_cli}, ${rep})"
    else
        line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${QUEUE_FILE}")
        if [ -z "${line}" ]; then
            echo "No queue entry for task ${SLURM_ARRAY_TASK_ID}, exiting"
            exit 0
        fi
        cond_cli="${line%%,*}"
        rep="${line##*,}"
    fi

    echo "============================================================"
    echo "Task ${SLURM_JOB_ID:-?}/${SLURM_ARRAY_TASK_ID}  ${cond_cli}  rep${rep}"
    echo "Node: $(hostname)"
    echo "GPU:  $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
    echo "Date: $(date)"
    echo "============================================================"

    # SIGUSR1 trap: SLURM sends this 90s before walltime (configured by
    # --signal=B:USR1@90 above). We sbatch a follow-up task that resumes
    # the SAME (cond, rep) and then exit gracefully so the partial blocks
    # already on disk are not corrupted by a forced kill.
    on_walltime() {
        echo "[$(date +%H:%M:%S)] SIGUSR1 received (walltime imminent)" \
             "— scheduling resume task for ${cond_cli} rep${rep}"
        # Kill the python child first so it stops writing.
        if [ -n "${PYPID:-}" ] && kill -0 "${PYPID}" 2>/dev/null; then
            kill -TERM "${PYPID}" 2>/dev/null || true
            wait "${PYPID}" 2>/dev/null || true
        fi
        sbatch \
            --array=0 \
            --export=ALL,RESUME_COND="${cond_cli}",RESUME_REP="${rep}" \
            "$0" \
            && echo "[$(date +%H:%M:%S)] resume task submitted." \
            || echo "[$(date +%H:%M:%S)] resume sbatch FAILED."
        exit 0
    }
    trap on_walltime SIGUSR1

    # Always pass --resume; if the rep dir is fresh, run_simulation.py
    # falls back to a clean grow_cubic init. --n-blocks caps the rep
    # length so it fits comfortably under the 48 h walltime.
    python scripts/run_simulation.py \
        --condition "${cond_cli}" \
        --replicate "${rep}" \
        --engine polychrom \
        --gpu 0 \
        --output "${RESULTS_BASE}" \
        --resume \
        --n-blocks "${BLOCKS_PER_REP}" \
        &
    PYPID=$!
    wait "${PYPID}"
    PY_RC=$?
    trap - SIGUSR1

    echo "Done: ${cond_cli} rep${rep} (python rc=${PY_RC}) — $(date)"
    exit "${PY_RC}"
fi

# ---------- Planner mode ----------
mkdir -p logs cluster
: > "${QUEUE_FILE}"

echo "============================================================"
echo "CATCH-UP PLANNER  ($(date))"
echo "Target frames/cond:    ${TARGET_FRAMES_PER_CONDITION}"
echo "Expected new rep size: ${EXPECTED_BLOCKS_PER_REP} blocks"
echo "                       ($((EXPECTED_BLOCKS_PER_REP * FRAMES_PER_BLOCK)) frames)"
echo "Include new conds:     ${INCLUDE_NEW_CONDITIONS}"
echo "============================================================"

declare -A current_frames
declare -A max_rep
declare -A n_reps
declare -A rep_summary
declare -A seen_reps

# Scan every rep dir (merged_X_repN or X_repN). Shard dirs (X_repN_shardM)
# do not end in `_repN` digits-only, so the regex below excludes them.
#
# Two passes so `merged_*` always wins over any bare-name leftover from
# pre-2026-04-27 (when merge_shards.py started prefixing). The bare dir
# is a stale partial of the same rep; counting it would under-report
# `current_frames` and cause MISSING_ONLY=0 runs to over-schedule.
shopt -s nullglob

scan_dir() {
    local d="$1"
    [ -d "$d" ] || return 0
    local name
    name=$(basename "$d")
    name="${name#merged_}"
    [[ "$name" =~ ^(.+)_rep([0-9]+)$ ]] || return 0
    local disk_name="${BASH_REMATCH[1]}"
    local rep="${BASH_REMATCH[2]}"

    [ -n "${DISK_TO_CONDITION[$disk_name]:-}" ] || return 0

    local key="${disk_name}|${rep}"
    if [ -n "${seen_reps[$key]:-}" ]; then
        # Only worth a warning if the skipped bare-name dir actually has
        # blocks (i.e. a real older copy, not an empty stub).
        local stale_blocks
        stale_blocks=$(ls "$d"/blocks_*.h5 2>/dev/null | wc -l)
        if [ "$stale_blocks" -gt 0 ]; then
            echo "WARN: bare-name leftover ${d} ($stale_blocks blocks) shadowed by merged_${disk_name}_rep${rep} — using merged_; bare dir is a candidate for cleanup" >&2
        fi
        return 0
    fi
    seen_reps[$key]=1

    local n_blocks n_frames
    n_blocks=$(ls "$d"/blocks_*.h5 2>/dev/null | wc -l)
    n_frames=$(( n_blocks * FRAMES_PER_BLOCK ))

    current_frames[$disk_name]=$(( ${current_frames[$disk_name]:-0} + n_frames ))
    local cur_max=${max_rep[$disk_name]:--1}
    if [ "$rep" -gt "$cur_max" ]; then max_rep[$disk_name]=$rep; fi
    n_reps[$disk_name]=$(( ${n_reps[$disk_name]:-0} + 1 ))
    rep_summary[$disk_name]+="    rep${rep}: ${n_blocks} blocks (${n_frames} frames)"$'\n'
}

# Pass 1: merged_* dirs (canonical post-refactor form)
for d in "${RESULTS_BASE}"/merged_*_rep[0-9]*; do scan_dir "$d"; done
# Pass 2: bare-name dirs (skipped if a merged_ counterpart was already seen)
for d in "${RESULTS_BASE}"/*_rep[0-9]*; do
    name=$(basename "$d")
    [[ "$name" == merged_* ]] && continue
    scan_dir "$d"
done

shopt -u nullglob

queue_size=0
for disk_name in $(printf '%s\n' "${!DISK_TO_CONDITION[@]}" | sort); do
    cond_cli="${DISK_TO_CONDITION[$disk_name]}"
    cur=${current_frames[$disk_name]:-0}

    if [ "$cur" -eq 0 ] && [ "${INCLUDE_NEW_CONDITIONS}" = "0" ]; then
        echo
        echo "[${disk_name}] (no data, skipped — INCLUDE_NEW_CONDITIONS=0)"
        continue
    fi
    if [ "$cur" -gt 0 ] && [ "${MISSING_ONLY}" = "1" ]; then
        echo
        echo "[${disk_name}] (already has data, skipped — MISSING_ONLY=1)"
        continue
    fi

    deficit=$(( TARGET_FRAMES_PER_CONDITION - cur ))
    [ "$deficit" -lt 0 ] && deficit=0
    new_frames=$(( EXPECTED_BLOCKS_PER_REP * FRAMES_PER_BLOCK ))
    n_new=0
    if [ "$deficit" -gt 0 ]; then
        n_new=$(( (deficit + new_frames - 1) / new_frames ))
    fi

    # Floor by MIN_REPS_PER_CONDITION (e.g. for MSD statistics): make sure
    # the condition ends up with at least that many reps total.
    existing_reps=${n_reps[$disk_name]:-0}
    min_extra=$(( MIN_REPS_PER_CONDITION - existing_reps ))
    [ "$min_extra" -lt 0 ] && min_extra=0
    if [ "$n_new" -lt "$min_extra" ]; then n_new=$min_extra; fi

    echo
    echo "[${disk_name}]  current=${cur} frames in ${existing_reps} reps,"
    echo "                deficit=${deficit}, min_extra_for_${MIN_REPS_PER_CONDITION}_reps=${min_extra}, new_reps=${n_new}"
    if [ -n "${rep_summary[$disk_name]:-}" ]; then
        echo -n "${rep_summary[$disk_name]}"
    fi

    existing_max=${max_rep[$disk_name]:--1}
    if [ "$existing_max" -lt 0 ]; then
        start_rep=1   # brand-new condition: no rep history to avoid, start at 1
    else
        start_rep=$(( existing_max + 1 ))
    fi

    for ((i = 0; i < n_new; i++)); do
        new_rep=$(( start_rep + i ))
        echo "${cond_cli},${new_rep}" >> "${QUEUE_FILE}"
        echo "    -> schedule: ${cond_cli} rep${new_rep}"
        queue_size=$(( queue_size + 1 ))
    done
done

echo
echo "============================================================"
echo "Queue size: ${queue_size}"
echo "Queue file: ${QUEUE_FILE}"
echo "============================================================"

if [ "${PLAN_ONLY:-0}" = "1" ]; then
    echo
    echo "PLAN_ONLY=1 — not submitting. Queue contents:"
    if [ -s "${QUEUE_FILE}" ]; then cat "${QUEUE_FILE}"; else echo "(empty)"; fi
    exit 0
fi
if [ "$queue_size" -eq 0 ]; then
    echo "Nothing to do. All conditions at target."
    exit 0
fi

echo
echo "Submitting SLURM array 0-$((queue_size - 1))%2 ..."
ARRAY_JID=$(sbatch --parsable --array=0-$((queue_size - 1))%2 "$0")
echo "  catchup array: ${ARRAY_JID}"

if [ "${CHAIN_FOLLOW_UPS}" != "1" ]; then
    echo
    echo "CHAIN_FOLLOW_UPS=0 — not chaining comparable_merge / analysis."
    echo "Run them manually after the array drains:"
    echo "  sbatch cluster/submit_comparable_merge.sh"
    echo "  sbatch cluster/submit_analysis_postmerge.sh"
    exit 0
fi

# Chain comparable_merge after the array — afterany so it fires even if a
# task hits walltime (the SIGUSR1 trap exits 0, but afterany doesn't care
# either way). Resume tasks spawned by the trap are separate jobs and
# typically don't fire because BLOCKS_PER_REP fits inside walltime; if they
# do, comparable_merge will warn for the affected reps and you can re-run
# it once the resume jobs drain.
COMP_JID=$(sbatch --parsable \
    --dependency=afterany:${ARRAY_JID} \
    cluster/submit_comparable_merge.sh)
echo "  comparable_merge: ${COMP_JID}  (after catchup array)"

# Then analysis on the comparable dirs (and refreshed regular pools).
ANA_JID=$(sbatch --parsable \
    --dependency=afterok:${COMP_JID} \
    cluster/submit_analysis_postmerge.sh)
echo "  analysis_postmerge: ${ANA_JID}  (after comparable_merge)"

echo
echo "Chain queued. Email notifications will fire on END/FAIL/TIME_LIMIT for"
echo "each step. Watch progress with:"
echo "  squeue -u \$USER -o '%10i %20j %2t %10M %20R'"
