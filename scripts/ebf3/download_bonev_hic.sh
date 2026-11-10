#!/bin/bash
#SBATCH --job-name=bonev_hic
#SBATCH --account=rc_dphansti_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=8
#SBATCH --mem=48g
#SBATCH --time=12:00:00
#SBATCH --output=data/logs/bonev_hic_%x_%a_%j.out
#SBATCH --array=0-2
# ============================================================================
# Bonev 2017 Hi-C valid pairs for one chromosome (mm10; CHROM env var, default chr7) -> .mcool
#
# Adapted from legacy/scripts/download_bonev_hic.sh (chr3 / Sox2). Changes:
#   * chr3 -> chr7; cell types ES, ncx_NPC, ncx_CN.
#   * Local files KEEP the remote `ncx_` prefix. GSE161259 reanalyses only the
#     neocortex-sorted NPC/CN (GSM2533835-42), NOT the in vitro NPC/CN that the
#     Bonev ChIP/RNA come from; the name must carry that distinction.
#   * The tabix / cload / zoomify commands are the legacy ones verbatim.
#
# Source: GEO GSE161259 (lldelisle HiCUP 0.6.1 reanalysis of Bonev 2017, mm10).
# Requires: tabix, cooler (both in the pipeline's cohesin_sim env).
#
# Usage (from the repo root):
#   CHROM=chr2 sbatch --export=ALL scripts/ebf3/download_bonev_hic.sh   # array: ES, ncx_NPC, ncx_CN
#   CHROM=chr2 bash scripts/ebf3/download_bonev_hic.sh ncx_CN     # one sample, interactively
# ============================================================================
set -euo pipefail

CHROM="${CHROM:-chr7}"   # export CHROM=chr2 for another chromosome
SAMPLES=(ES ncx_NPC ncx_CN)
if [ $# -ge 1 ]; then SAMPLE="$1"; else SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}"; fi
case "${SAMPLE}" in ES|ncx_NPC|ncx_CN) ;; *) echo "sample must be ES, ncx_NPC or ncx_CN" >&2; exit 1 ;; esac

SIMENV="${SIMENV:-.snakemake/conda/0544ae75e3043e68c45ece7022cee694_}"
export PATH="${SIMENV}/bin:${PATH}"

DATA_DIR="data/mcool"
PAIRS_DIR="data/raw/gse161259"          # untrusted downloads in their own folder
mkdir -p "${DATA_DIR}" "${PAIRS_DIR}" data/logs

FTP_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161259/suppl"
BIN=1
NPROC="${NPROC:-${SLURM_CPUS_PER_TASK:-8}}"
MAX_SPLIT="${MAX_SPLIT:-${NPROC}}"

# UCSC mm10 chrom sizes (data/raw/ucsc/mm10.chrom.sizes), primary chromosomes only,
# in the same order the legacy script used.
SIZES_FILE="data/mm10.chrom.sizes"
if [ ! -s "${SIZES_FILE}" ]; then
    grep -E "^chr([0-9]+|X|Y)\s" data/raw/ucsc/mm10.chrom.sizes | sort -V > "${SIZES_FILE}"
fi

REMOTE_NAME="GSE161259_${SAMPLE}_${CHROM}.validPairs.csort.txt.gz"
PAIRS_FILE="${PAIRS_DIR}/${REMOTE_NAME}"
if [ ! -s "${PAIRS_FILE}" ]; then
    echo "[${SAMPLE}] downloading ${REMOTE_NAME}"
    wget -q -O "${PAIRS_FILE}.part" "${FTP_BASE}/${REMOTE_NAME}"
    mv "${PAIRS_FILE}.part" "${PAIRS_FILE}"
fi
gzip -t "${PAIRS_FILE}"

if [ ! -f "${PAIRS_FILE}.tbi" ]; then
    echo "[${SAMPLE}] tabix"
    tabix -s 3 -b 4 -e 4 "${PAIRS_FILE}"
fi

RAW_FILE="${DATA_DIR}/${SAMPLE}_${CHROM}_raw.${BIN}kb.cool"
MCOOL_FILE="${DATA_DIR}/${SAMPLE}_${CHROM}.mcool"
if [ -f "${RAW_FILE}" ] && ! cooler info "${RAW_FILE}" 2>/dev/null | grep -q "bin-size"; then
    rm -f "${RAW_FILE}"
fi
rm -f "${MCOOL_FILE}"

if [ ! -f "${RAW_FILE}" ]; then
    echo "[${SAMPLE}] cooler cload tabix (nproc=${NPROC}, max-split=${MAX_SPLIT})"
    cooler cload tabix \
        --nproc "${NPROC}" \
        -c2 7 -p2 8 -s "${MAX_SPLIT}" \
        --assembly mm10 \
        "${SIZES_FILE}:${BIN}000" "${PAIRS_FILE}" "${RAW_FILE}"
fi

echo "[${SAMPLE}] cooler zoomify + balance"
cooler zoomify \
    --nproc "${NPROC}" \
    --balance \
    --balance-args "--cis-only" \
    --resolutions 1000,5000,10000,25000,50000,100000 \
    --out "${MCOOL_FILE}" \
    "${RAW_FILE}"

cooler ls "${MCOOL_FILE}"
echo "[${SAMPLE}] done: ${MCOOL_FILE}"
