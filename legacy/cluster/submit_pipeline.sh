#!/bin/bash
# ============================================================================
# cluster/submit_pipeline.sh - END-TO-END PIPELINE LAUNCHER
#
# Submits the full simulation -> merge -> analysis chain as ONE dependency
# tree. The analysis stage is itself a CHAIN of resume-capable jobs linked by
# --dependency=afterany, so wall-time hits, OOMs, or transient cluster
# faults don't break the pipeline - the next link wakes up automatically and
# picks up where the previous one stopped (--reuse-heavy + --skip-existing).
#
# This is the recommended entry point for a fresh sweep. Run once, walk away.
#
# Stages
# ------
#   1. cluster/submit_polychrom_multigpu.sh        Tier 2 GPU sims (8 array tasks)
#   2. cluster/submit_merge_shards.sh              afterok step 1 -> merged_*/
#   3. cluster/submit_analysis_postmerge.sh        afterok step 2  (link 1 of N)
#   3b. ...same script, RESUME=1, SKIP_EXISTING=1  afterany previous link
#       (repeated --resume-links N times; each link short-circuits if the
#        previous one already finished, so extra links are essentially free)
#
# Usage
# -----
#   bash cluster/submit_pipeline.sh                     # full chain, 3 analysis links
#   bash cluster/submit_pipeline.sh --skip-sim          # skip sim, start at merge
#   bash cluster/submit_pipeline.sh --skip-merge        # skip merge, start at analysis
#   bash cluster/submit_pipeline.sh --resume-links 5    # 5 analysis chain links
#   bash cluster/submit_pipeline.sh --analysis-script cluster/submit_analysis_all_16cpu.sh
#
#   # Attach the chain to an already-running job (e.g. you submitted sim
#   # manually and want this script to pick up at merge after it ends):
#   bash cluster/submit_pipeline.sh --skip-sim --after 9457828
#
# Stop the chain at any point:
#   scancel <jobid>          # cancel one link; afterany still fires the next one
#   scancel -u $USER         # cancel everything you own
#
# Notes
# -----
# - The first analysis link runs cold (no SKIP_EXISTING). Subsequent links
#   set RESUME=1 and SKIP_EXISTING=1 so they short-circuit anything the
#   previous link already produced. Net: extra links cost seconds, not hours.
# - --dependency=afterany is intentional. afterok would stop the chain when
#   a link times out (which is exactly when we want the next link to fire).
#   Use scancel to actually stop the chain.
# - All env-var toggles for the underlying scripts (CLEANUP, SKIP_MSD,
#   SKIP_APA, etc.) are passed through unchanged via the user environment.
# ============================================================================

set -euo pipefail

# Defaults
N_RESUME_LINKS=3
ANALYSIS_SCRIPT="cluster/submit_analysis_postmerge.sh"
SIM_SCRIPT="cluster/submit_polychrom_multigpu.sh"
MERGE_SCRIPT="cluster/submit_merge_shards.sh"
SKIP_SIM=0
SKIP_MERGE=0
ATTACH_TO=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --skip-sim)        SKIP_SIM=1; shift ;;
        --skip-merge)      SKIP_MERGE=1; shift ;;
        --resume-links)    N_RESUME_LINKS="$2"; shift 2 ;;
        --analysis-script) ANALYSIS_SCRIPT="$2"; shift 2 ;;
        --sim-script)      SIM_SCRIPT="$2"; shift 2 ;;
        --merge-script)    MERGE_SCRIPT="$2"; shift 2 ;;
        --after)           ATTACH_TO="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,55p' "$0"; exit 0 ;;
        *)
            echo "Unknown arg: $1" >&2
            echo "Run with --help for usage." >&2
            exit 1 ;;
    esac
done

if ! [[ "${N_RESUME_LINKS}" =~ ^[0-9]+$ ]] || [ "${N_RESUME_LINKS}" -lt 1 ]; then
    echo "ERROR: --resume-links must be a positive integer (got '${N_RESUME_LINKS}')." >&2
    exit 1
fi

cd "$(dirname "$0")/.."
mkdir -p logs

for s in "${ANALYSIS_SCRIPT}" "${SIM_SCRIPT}" "${MERGE_SCRIPT}"; do
    [ -f "$s" ] || { echo "ERROR: missing wrapper $s" >&2; exit 2; }
done

echo "================================================================"
echo " PIPELINE SUBMIT - $(date)"
echo "================================================================"
echo "  sim script:        ${SIM_SCRIPT}     skip=${SKIP_SIM}"
echo "  merge script:      ${MERGE_SCRIPT}   skip=${SKIP_MERGE}"
echo "  analysis script:   ${ANALYSIS_SCRIPT}"
echo "  resume links:      ${N_RESUME_LINKS}"
echo "  attach to:         ${ATTACH_TO:-<none>}"
echo "================================================================"

PREV_JID="${ATTACH_TO}"

# --- Stage 1: simulation ----------------------------------------------------
if [ "${SKIP_SIM}" != "1" ]; then
    if [ -n "${PREV_JID}" ]; then
        SIM_JID=$(sbatch --parsable --dependency=afterok:"${PREV_JID}" "${SIM_SCRIPT}")
    else
        SIM_JID=$(sbatch --parsable "${SIM_SCRIPT}")
    fi
    echo "  [1] sim         JID=${SIM_JID}   (${SIM_SCRIPT})"
    PREV_JID="${SIM_JID}"
else
    echo "  [1] sim         SKIPPED"
fi

# --- Stage 2: merge ---------------------------------------------------------
if [ "${SKIP_MERGE}" != "1" ]; then
    if [ -n "${PREV_JID}" ]; then
        MERGE_JID=$(sbatch --parsable --dependency=afterok:"${PREV_JID}" "${MERGE_SCRIPT}")
    else
        MERGE_JID=$(sbatch --parsable "${MERGE_SCRIPT}")
    fi
    echo "  [2] merge       JID=${MERGE_JID}   (${MERGE_SCRIPT})"
    PREV_JID="${MERGE_JID}"
else
    echo "  [2] merge       SKIPPED"
fi

# --- Stage 3: analysis chain (N links) -------------------------------------
ANALYSIS_JIDS=()
for i in $(seq 1 "${N_RESUME_LINKS}"); do
    if [ "$i" -eq 1 ]; then
        # First link: cold run. afterok on the merge stage if that ran.
        if [ -n "${PREV_JID}" ]; then
            JID=$(sbatch --parsable --dependency=afterok:"${PREV_JID}" "${ANALYSIS_SCRIPT}")
        else
            JID=$(sbatch --parsable "${ANALYSIS_SCRIPT}")
        fi
    else
        # Subsequent link: afterany on previous, RESUME=1 + SKIP_EXISTING=1.
        # We pass these via --export so the wrapper picks them up.
        JID=$(sbatch --parsable --dependency=afterany:"${PREV_JID}" \
              --export=ALL,RESUME=1,SKIP_EXISTING=1 \
              "${ANALYSIS_SCRIPT}")
    fi
    ANALYSIS_JIDS+=("${JID}")
    echo "  [3.${i}] analysis JID=${JID}   $( [ $i -eq 1 ] && echo "(cold)" || echo "(resume, afterany ${PREV_JID})" )"
    PREV_JID="${JID}"
done

echo "================================================================"
echo "All submitted. Final JID = ${PREV_JID}"
echo
echo "Useful commands:"
echo "  squeue -u \$USER --format='%.10i %.20j %.10T %.10M %.20R'"
echo "  scontrol show job ${ANALYSIS_JIDS[0]}"
echo "  tail -f logs/analysis_postmerge_${ANALYSIS_JIDS[0]}.out"
echo "================================================================"
