#!/bin/bash
# ============================================================================
# check_gpu_util.sh — live GPU utilisation snapshot for currently-running
# catchup tasks. Use this to decide whether asking for a higher 2-GPU cap
# (or enabling MPS) is worth it.
#
# How to read the output:
#   utilization.gpu  — % of the time the GPU was busy in the last sample.
#                      80-95% sustained = GPU well-fed; gain from MPS or
#                      more concurrent processes is small (~5-15%).
#                      30-60% = GPU underutilised; MPS / more processes
#                      could ~2x throughput.
#   memory.used      — VRAM in use. polychrom on a 70k-monomer system
#                      typically sits at 1-3 GB / 16-32 GB on V100,
#                      so fitting 2-3 sims on one GPU is feasible
#                      memory-wise (the question is just compute).
#
# Usage
# -----
#   bash cluster/check_gpu_util.sh                 # auto-find the catchup node
#   bash cluster/check_gpu_util.sh dgx01           # explicit node
#   POLL=5 bash cluster/check_gpu_util.sh          # refresh every 5s (default 10)
#   N_POLLS=6 bash cluster/check_gpu_util.sh       # collect 6 snapshots and exit
#
# Run from a cluster login node (it ssh's into the compute node).
# ============================================================================

set -euo pipefail

POLL="${POLL:-10}"
N_POLLS="${N_POLLS:-0}"   # 0 = run forever (Ctrl-C to stop)

if [ $# -ge 1 ]; then
    node="$1"
else
    node=$(squeue -u "$USER" -h -t R --name=catchup -o '%R' 2>/dev/null \
           | grep -v '^(' | sort -u | head -1 | tr -d '[:space:]')
    if [ -z "${node}" ]; then
        echo "No running catchup task found. Current state:"
        squeue -u "$USER" --name=catchup -o '%10i %2t %20R' 2>/dev/null \
            | head -20
        echo
        echo "Re-run once a task is in state R, or pass a node explicitly:"
        echo "  bash cluster/check_gpu_util.sh dgx01"
        exit 1
    fi
fi

echo "Polling nvidia-smi on ${node} every ${POLL}s "
[ "${N_POLLS}" -gt 0 ] && echo "(${N_POLLS} samples then exit)" || echo "(Ctrl-C to stop)"
echo "----------------------------------------------------------------"
printf '%-9s %-6s %-9s %-12s %-5s\n' \
    "time" "gpu" "util%" "vram_used_MB" "tempC"
echo "----------------------------------------------------------------"

count=0
while :; do
    snapshot=$(ssh -o BatchMode=yes -o StrictHostKeyChecking=no "${node}" \
        "nvidia-smi --query-gpu=index,utilization.gpu,memory.used,temperature.gpu \
                    --format=csv,noheader,nounits" 2>/dev/null) || {
        echo "ssh ${node} failed; aborting."
        exit 1
    }
    ts=$(date '+%H:%M:%S')
    while IFS=, read -r idx util mem temp; do
        printf '%-9s %-6s %-9s %-12s %-5s\n' \
            "$ts" \
            "$(echo "$idx" | tr -d '[:space:]')" \
            "$(echo "$util" | tr -d '[:space:]')" \
            "$(echo "$mem" | tr -d '[:space:]')" \
            "$(echo "$temp" | tr -d '[:space:]')"
    done <<< "$snapshot"
    count=$((count + 1))
    if [ "${N_POLLS}" -gt 0 ] && [ "${count}" -ge "${N_POLLS}" ]; then break; fi
    sleep "${POLL}"
done
