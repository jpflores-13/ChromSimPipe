#!/bin/bash
# ============================================================================
# CLUSTER PARTITION & RESOURCE DISCOVERY
#
# Run this interactively on the cluster (NOT via sbatch) to find the maximum
# CPUs and memory you can request, and which partitions are available.
#
# Usage:
#   bash cluster/check_partitions.sh
# ============================================================================

SEP="============================================================"

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "1. PARTITION SUMMARY  (name / state / max-time / total-nodes / CPUs)"
echo "$SEP"
# Note: some SLURM versions don't support printf-width specifiers (%-20P).
# Using plain %X codes here for compatibility.
sinfo -o "%P %a %l %D %C %m" 2>/dev/null \
    | column -t \
    || { echo "  (sinfo not available — run this on the cluster)"; exit 1; }
echo ""
echo "  CPU column format: allocated/idle/other/total"

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "2. NODES: actual hardware per node (cores / RAM / GPUs)"
echo "$SEP"
sinfo -N -o "%N %P %t %c %m %G" 2>/dev/null | column -t | sort -k2,2
echo ""
echo "  Memory shown in MB.  Divide by 1024 for GB."

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "3. MAXIMUM CPUs AND MEMORY YOU CAN REQUEST PER JOB"
echo "   (per-node hardware ceiling, before any QOS limits)"
echo "$SEP"
echo ""
echo "  By partition:"
# For each partition, find the node with the most CPUs and most RAM
sinfo -N -o "%P %c %m" 2>/dev/null \
    | tail -n +2 \
    | awk '
    {
        part=$1; cpus=$2; mem=$3
        if (cpus > max_cpu[part]) max_cpu[part]=cpus
        if (mem  > max_mem[part]) max_mem[part]=mem
    }
    END {
        printf "  %-25s %10s %15s\n", "Partition", "Max CPUs", "Max Memory (GB)"
        printf "  %-25s %10s %15s\n", "---------", "--------", "---------------"
        for (p in max_cpu)
            printf "  %-25s %10d %15.0f\n", p, max_cpu[p], max_mem[p]/1024
    }' \
    | sort -k1

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "4. YOUR QOS LIMITS  (hard caps set by the cluster admin for your account)"
echo "$SEP"
echo ""
echo "  QOS table (MaxCPUsPU = max CPUs per user across all jobs):"
sacctmgr show qos \
    format=Name%-20,MaxCPUsPU,MaxGRESPU%-20,MaxJobsPU,MaxSubmitPU,MaxWallDurationPerJob \
    2>/dev/null || echo "  (sacctmgr not available)"

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "5. YOUR ACCOUNT ASSOCIATIONS  (partitions + QOS you are allowed to use)"
echo "$SEP"
sacctmgr show association user=$USER \
    format=Account%-20,Partition%-20,QOS%-30,MaxCPUs,MaxJobs,MaxSubmit \
    2>/dev/null || echo "  (sacctmgr not available)"

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "6. CURRENT USAGE  (how many CPUs/GPUs you are using right now)"
echo "$SEP"
squeue -u $USER \
    -o "%.10i %-20j %-10T %-8C %-12m %-12l %P %R" \
    2>/dev/null || echo "  (no jobs or squeue unavailable)"
echo ""
USED_CPUS=$(squeue -u $USER -h -o "%C" 2>/dev/null | awk '{s+=$1} END{print s+0}')
echo "  Total CPUs currently in use by you: ${USED_CPUS}"

# ---------------------------------------------------------------------------
echo ""
echo "$SEP"
echo "7. QUICK RECOMMENDATION FOR submit_analysis_cpu.sh"
echo "$SEP"
echo ""
echo "  1. Find the largest CPU-only partition in section 1 (no GPU column)."
echo "  2. From section 3, read its Max CPUs and Max Memory."
echo "  3. Check section 4 for any per-user QOS cap (MaxCPUsPU) — take the"
echo "     LOWER of the hardware max and the QOS cap."
echo "  4. Edit submit_analysis_cpu.sh:"
echo ""
echo "       #SBATCH --partition=<partition_name>"
echo "       #SBATCH --cpus-per-task=<max_cpus>"
echo "       #SBATCH --mem=<floor(max_cpus/4) * 7>G   # 7 GB per concurrent dir"
echo ""
echo "  Example: 80 CPUs, no QOS cap → --cpus-per-task=80 --mem=150G"
echo "  Example: 32 CPUs, QOS cap 24 → --cpus-per-task=24 --mem=44G"
echo ""
echo "$SEP"
