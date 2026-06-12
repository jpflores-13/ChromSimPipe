#!/bin/bash
# Submit a CHAIN of analysis sbatch jobs, each depending on the previous
# via --dependency=afterany so the next link fires regardless of how the
# previous one ended (COMPLETED / TIMEOUT / FAILED / CANCELLED).
#
# Each link runs with RESUME=1 + SKIP_EXISTING=1, so it skips merged
# directories whose {condition}_*blk_rep*_contact_map.npy already exists
# in results/analysis/, and reuses cached contact maps + MSD outputs for
# the rest. Net effect: a single chain of K jobs covers K x walltime of
# wall-clock progress on the same sweep.
#
# This is the standard answer for sweeps that don't finish in 48 h.
#
# Usage
# -----
#   # Chain of 3 jobs starting fresh (recommended starting point)
#   bash cluster/submit_analysis_resume_chain.sh
#
#   # Chain of N=5 jobs
#   bash cluster/submit_analysis_resume_chain.sh 5
#
#   # Chain of 3 jobs that ATTACHES to an already-running job
#   # (the new chain only starts after job 9457828 ends, in any state).
#   bash cluster/submit_analysis_resume_chain.sh 3 9457828
#
#   # Use the 16-CPU sequential variant instead of the 32-CPU 2x parallel
#   bash cluster/submit_analysis_resume_chain.sh 5 "" cluster/submit_analysis_all_16cpu.sh
#
# Stop the chain at any point with:
#   scancel <jobid>      # cancel one link; subsequent links still fire on afterany
#   scancel -u $USER     # cancel everything you own
#
# To prevent later links from firing when an earlier link has clearly
# finished the work, the chain is fire-and-forget: every link will run.
# That is intentional - subsequent links have nothing to do (skip-existing
# + reuse-heavy short-circuits everything) and exit in seconds, so the
# extra runs are essentially free.

set -euo pipefail

N=${1:-3}
INITIAL_DEP=${2:-}
WRAPPER=${3:-cluster/submit_analysis_all_32cpu.sh}

if ! [[ "$N" =~ ^[0-9]+$ ]] || [ "$N" -lt 1 ] || [ "$N" -gt 20 ]; then
    echo "ERROR: chain length must be an integer in [1, 20]; got '${N}'" >&2
    exit 1
fi
if [ ! -f "$WRAPPER" ]; then
    echo "ERROR: wrapper not found: ${WRAPPER}" >&2
    exit 1
fi

mkdir -p logs

EXPORT_ARGS="--export=ALL,RESUME=1,SKIP_EXISTING=1"
prev_dep=""
if [ -n "$INITIAL_DEP" ]; then
    prev_dep="--dependency=afterany:${INITIAL_DEP}"
    echo "Chaining ${N} links after existing job ${INITIAL_DEP}."
else
    echo "Submitting fresh chain of ${N} links."
fi

declare -a chain
for i in $(seq 1 "$N"); do
    # First link: optional initial dependency; subsequent links: depend on previous link.
    if [ "$i" -eq 1 ]; then
        deps="${prev_dep}"
    else
        deps="--dependency=afterany:${chain[$((i-2))]}"
    fi

    jid=$(sbatch --parsable ${deps} ${EXPORT_ARGS} "$WRAPPER")
    chain+=("$jid")
    echo "  link ${i}/${N}: ${jid}    deps='${deps}'"
done

echo
echo "Chain submitted. ${N} job IDs: ${chain[*]}"
echo
echo "Watch the queue:"
echo "  squeue -u \$USER -o '%10i %20j %2t %10M %20R'"
echo
echo "Cancel the rest of the chain (keeping the running link):"
echo "  scancel ${chain[*]:1}"
