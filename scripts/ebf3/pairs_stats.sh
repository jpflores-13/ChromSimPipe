#!/bin/bash
#SBATCH --job-name=pairs_stats
#SBATCH --account=rc_dphansti_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=4g
#SBATCH --time=3:00:00
#SBATCH --array=0-2
#SBATCH --output=data/logs/pairs_stats_%a_%j.out
## Ebf3 Stage 1d: valid-pair counts in the GSE161259 chr7 files (cis = both ends chr7).
set -euo pipefail
S=(ES ncx_NPC ncx_CN); s=${S[$SLURM_ARRAY_TASK_ID]}
zcat data/raw/gse161259/GSE161259_${s}_chr7.validPairs.csort.txt.gz | awk -v s=$s '
  {n++; if ($3==$7) {c++; d=$8-$4; if (d<0) d=-d; if (d<1000) l1++; else if (d<20000) l20++; else g20++} else t++}
  END {printf "{\"valid_pairs_in_file\": %d, \"cis_chr7\": %d, \"trans\": %d, \"cis_lt_1kb\": %d, \"cis_1_20kb\": %d, \"cis_ge_20kb\": %d}\n", n, c, t, l1, l20, g20}' > results/ebf3/pairs_stats_${s}.json
cat results/ebf3/pairs_stats_${s}.json
