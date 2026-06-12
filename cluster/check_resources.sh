#!/bin/bash
# ============================================================================
# check_resources.sh — comprehensive resource-usage snapshot for the running
# catchup tasks. Run on dgx01 (or wherever your tasks are scheduled).
#
# Captures:
#   1. Per-GPU utilization sustained over 60 s (24 samples × 2.5 s)
#   2. VRAM and PID per CUDA process
#   3. CPU%, RAM (RSS), elapsed time per *user* python process
#   4. Map of (GPU index → user) for at-a-glance ownership
#
# Run it WHILE catchup tasks are mid-flight (i.e. during steady-state
# polymer dynamics, not during the brief LEF pre-compute that opens each
# rep). 60 s of GPU sampling is enough to smooth nvidia-smi's instantaneous
# noise while staying short.
#
# Usage
# -----
#   # Already on dgx01:
#   bash cluster/check_resources.sh
#
#   # From the login node:
#   ssh dgx01 bash /beegfs/scratch/ric.gabriele/ric.gabriele/runs/cohesin_sim/cluster/check_resources.sh
#
#   # Tweak the GPU sampling window (default 60 s = 24 × 2.5 s):
#   GPU_SAMPLES=12 GPU_PERIOD=2 bash cluster/check_resources.sh
# ============================================================================

set -euo pipefail

GPU_SAMPLES="${GPU_SAMPLES:-24}"
GPU_PERIOD="${GPU_PERIOD:-2}"   # nvidia-smi dmon takes integer seconds (-d)

echo "============================================================"
echo "RESOURCE SNAPSHOT  ($(date))"
echo "host: $(hostname)   user: $USER"
echo "============================================================"
echo

# --- 1. GPU utilization over a window ---
echo "[1/4] Per-GPU utilization, ${GPU_SAMPLES} samples × ${GPU_PERIOD}s"
echo "      (average across the window is the steady-state number to trust)"
echo
nvidia-smi dmon -s u -c "${GPU_SAMPLES}" -d "${GPU_PERIOD}" 2>/dev/null \
    | awk 'NR<=2 {print; next}
           NR>2 {
               for (i=2; i<=NF; i++) sum[i]+=$i; cnt[i]++; print
           }
           END {
               printf "\navg utility per GPU index across %d samples:\n", NR-2
               for (i=2; i<=length(sum)+1; i++) {
                   if (sum[i] != "") printf "  gpu %d: ~%.1f %% (sm)\n", i-2, sum[i]/(NR-2)
               }
           }' \
    || echo "(nvidia-smi dmon failed)"

echo
echo "------------------------------------------------------------"

# --- 2. VRAM + PID per CUDA process ---
echo "[2/4] CUDA processes (per-process VRAM)"
echo
nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name,used_memory \
           --format=csv,noheader,nounits 2>/dev/null \
| while IFS=, read -r uuid pid name mem; do
    pid=${pid// /}; uuid=${uuid// /}; name=${name// /}; mem=${mem// /};
    user=$(ps -o user= -p "$pid" 2>/dev/null | tr -d '[:space:]')
    gpu=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader \
          | awk -v u="$uuid" -F',' '{gsub(/ /,""); if($2==u) print $1}')
    [ "$user" = "$USER" ] && mark="  <-- you" || mark=""
    printf "  GPU %-2s  PID %-8s  user %-12s  %-15s  %s MiB%s\n" \
        "$gpu" "$pid" "$user" "$name" "$mem" "$mark"
done

echo
echo "------------------------------------------------------------"

# --- 3. CPU% + RSS for the user's python processes ---
echo "[3/4] Your python processes (CPU%, RSS, elapsed)"
echo
ps -u "$USER" -o pid,pcpu,pmem,rss,etime,cmd 2>/dev/null \
    | awk 'NR==1 || /python.*run_simulation/'

echo
echo "------------------------------------------------------------"

# --- 4. Quick summary ---
echo "[4/4] Allocated-vs-used summary"
echo
my_python_count=$(pgrep -u "$USER" -f 'run_simulation.py' | wc -l)
my_gpu_count=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader \
               | tr -d ' ' \
               | while read pid; do
                   ps -o user= -p "$pid" 2>/dev/null | tr -d '[:space:]'
               done | grep -c "^${USER}$" || true)
total_rss_kb=$(ps -u "$USER" -o rss= -C python 2>/dev/null | awk '{s+=$1} END{print s+0}')
total_rss_gb=$(awk -v r="$total_rss_kb" 'BEGIN{printf "%.2f", r/1024/1024}')

echo "  python processes you own:        ${my_python_count}"
echo "  CUDA contexts you hold on dgx01: ${my_gpu_count}"
echo "  total Python RSS (RAM, GiB):     ${total_rss_gb}"
echo
echo "Decision hints for SIMS_PER_TASK (in-task multi-sim):"
echo "  - Sustained per-GPU sm% from [1] ÷ # of your sims = per-sim load"
echo "  - If per-sim GPU load <= 30%, expect ~2x speedup with SIMS_PER_TASK=2"
echo "  - If per-sim GPU load <= 15%, expect ~3x with SIMS_PER_TASK=4"
echo "  - SLURM allocates --cpus-per-task=4 and --mem=16G per task in"
echo "    submit_catch_up.sh; check that total RSS / N_sims fits in those."
