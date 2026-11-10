#!/bin/bash
#SBATCH --job-name=hic_halves
#SBATCH --account=rc_dphansti_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=48g
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=data/logs/hic_halves_%a_%j.out
# ============================================================================
# Ebf3 Stage 6 replicate ceiling: split each GSE161259 chr7 valid-pairs file
# into two random halves (each pair assigned by a seeded coin flip), and build
# one .mcool per half with the SAME cload / zoomify settings as Stage 1a.
# Every metric computed between halves is the ceiling no model may exceed.
#   sbatch scripts/ebf3/split_half_coolers.sh
# ============================================================================
set -euo pipefail
S=(ES ncx_NPC ncx_CN); s=${S[${SLURM_ARRAY_TASK_ID:-0}]}
SIMENV="${SIMENV:-.snakemake/conda/0544ae75e3043e68c45ece7022cee694_}"
export PATH="${SIMENV}/bin:${PATH}"
NPROC="${SLURM_CPUS_PER_TASK:-8}"
IN="data/raw/gse161259/GSE161259_${s}_chr7.validPairs.csort.txt.gz"
OUT="data/mcool/halves"; mkdir -p "${OUT}"
SIZES="data/mm10.chrom.sizes"

# Coin flip with a fixed seed per sample; input is position-sorted, so each
# half stays sorted and can be bgzipped + tabix-indexed like the full file.
zcat "${IN}" | awk -v seed=20261110 -v a="${OUT}/${s}_A.pairs" -v b="${OUT}/${s}_B.pairs" \
    'BEGIN{srand(seed)} {if (rand() < 0.5) print > a; else print > b}'
for h in A B; do
    bgzip -f -@ "${NPROC}" "${OUT}/${s}_${h}.pairs"
    tabix -f -s 3 -b 4 -e 4 "${OUT}/${s}_${h}.pairs.gz"
    cooler cload tabix --nproc "${NPROC}" -c2 7 -p2 8 -s "${NPROC}" --assembly mm10 \
        "${SIZES}:1000" "${OUT}/${s}_${h}.pairs.gz" "${OUT}/${s}_${h}_raw.1kb.cool"
    cooler zoomify --nproc "${NPROC}" --balance --balance-args "--cis-only" \
        --resolutions 1000,5000,10000,25000,50000,100000 \
        --out "${OUT}/${s}_${h}_chr7.mcool" "${OUT}/${s}_${h}_raw.1kb.cool"
    rm -f "${OUT}/${s}_${h}_raw.1kb.cool"
done
echo "[${s}] halves: $(for h in A B; do echo -n "${h}=$(cooler info ${OUT}/${s}_${h}_chr7.mcool::resolutions/1000 | grep '"sum"') "; done)"
