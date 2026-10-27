#!/bin/bash
# ============================================================================
# Download Bonev 2017 Hi-C valid pairs (mm10) from GEO GSE161259
# and convert to cool/mcool for the Sox2 locus (chr3:34-36 Mb).
#
# Source: lldelisle/Hi-C_reanalysis_Bonev_2017
# Only downloads chr3 (the Sox2 chromosome) to save time and space.
#
# Requirements: wget, tabix (htslib), cooler
#   conda install -c bioconda cooler htslib
#
# Usage (local, serial over both cell types):
#   bash scripts/download_bonev_hic.sh
#
# Usage (local, one cell type only — run two shells in parallel):
#   bash scripts/download_bonev_hic.sh ES &
#   bash scripts/download_bonev_hic.sh CN &
#   wait
#
# Usage (cluster, parallel ES + CN jobs):
#   sbatch --export=CELLTYPE=ES scripts/sbatch_bonev_hic.sbatch
#   sbatch --export=CELLTYPE=CN scripts/sbatch_bonev_hic.sbatch
#
# Knobs (env vars):
#   NPROC   number of workers for cooler cload/zoomify (default: SLURM
#           cpus-per-task, else 8).  More workers ≈ proportionally faster
#           binning and balancing for chr3.
# ============================================================================

set -euo pipefail

DATA_DIR="data/mcool"
PAIRS_DIR="data/valid_pairs"
mkdir -p "${DATA_DIR}" "${PAIRS_DIR}"

FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161259/suppl"

# Resolution for cool files (1 kb to match simulation)
BIN=1

# --- Which cell type(s) to process ---------------------------------------
# First positional arg selects a single cell type; no arg → both (serial).
if [ $# -ge 1 ]; then
    case "$1" in
        ES|CN) CELLTYPES=("$1") ;;
        *) echo "ERROR: cell type must be ES or CN (got '$1')" >&2; exit 1 ;;
    esac
else
    CELLTYPES=(ES CN)
fi

# --- Parallelism ---------------------------------------------------------
# Under SLURM we pick up the allocated cpus-per-task automatically; outside
# SLURM we default to 8.  MAX_SPLIT controls how many chunks chr3 is cut
# into for `cooler cload tabix` — with only 2 chunks, --nproc > 2 has
# nothing to do.  We tie it to NPROC so workers stay busy.
NPROC="${NPROC:-${SLURM_CPUS_PER_TASK:-8}}"
MAX_SPLIT="${MAX_SPLIT:-${NPROC}}"
echo "Parallelism: NPROC=${NPROC}  MAX_SPLIT=${MAX_SPLIT}  cell types: ${CELLTYPES[*]}"

# mm10 chromosome sizes (chr3 only needed, but cooler needs all)
SIZES_FILE="data/mm10.chrom.sizes"
if [ ! -f "${SIZES_FILE}" ]; then
    echo "Creating mm10 chrom sizes file..."
    cat > "${SIZES_FILE}" << 'SIZES'
chr1	195471971
chr2	182113224
chr3	160039680
chr4	156508116
chr5	151834684
chr6	149736546
chr7	145441459
chr8	129401213
chr9	124595110
chr10	130694993
chr11	122082543
chr12	120129022
chr13	120421639
chr14	124902244
chr15	104043685
chr16	98207768
chr17	94987271
chr18	90702639
chr19	61431566
chrX	171031299
chrY	91744698
SIZES
fi

# --- Download valid pairs for chr3 ---
for CELLTYPE in "${CELLTYPES[@]}"; do
    PAIRS_FILE="${PAIRS_DIR}/GSE161259_${CELLTYPE}_chr3.validPairs.csort.txt.gz"

    # Remote filename quirk on the lldelisle re-upload: ES is bare, but the
    # neural samples (CN, NPC) have an `ncx_` (neocortex tissue) prefix.
    # We always save with the clean "<CELLTYPE>_chr3" name locally so that
    # downstream steps don't have to know about the upload quirk.
    case "${CELLTYPE}" in
        ES)        REMOTE_NAME="GSE161259_${CELLTYPE}_chr3.validPairs.csort.txt.gz" ;;
        CN|NPC)    REMOTE_NAME="GSE161259_ncx_${CELLTYPE}_chr3.validPairs.csort.txt.gz" ;;
        *)         REMOTE_NAME="GSE161259_${CELLTYPE}_chr3.validPairs.csort.txt.gz" ;;
    esac

    if [ -f "${PAIRS_FILE}" ]; then
        echo "[${CELLTYPE}] Valid pairs already downloaded: ${PAIRS_FILE}"
    else
        echo "[${CELLTYPE}] Downloading chr3 valid pairs (remote: ${REMOTE_NAME})..."
        # Stream to .part and rename on success so a 404 / interrupted
        # download doesn't leave a truncated file that later runs would
        # mistake for a complete one.
        wget -O "${PAIRS_FILE}.part" "${FTP_BASE}/${REMOTE_NAME}"
        mv "${PAIRS_FILE}.part" "${PAIRS_FILE}"
    fi

    # --- Index with tabix ---
    if [ ! -f "${PAIRS_FILE}.tbi" ]; then
        echo "[${CELLTYPE}] Indexing with tabix..."
        tabix -s 3 -b 4 -e 4 "${PAIRS_FILE}"
    fi

    # --- Safety: remove any partial / incomplete outputs from earlier ---
    # A cload/balance/zoomify crash can leave a header-only HDF5 file
    # behind that the `if [ ! -f ... ]` guards would treat as "already
    # done".  We prefer idempotent re-runs: only the raw cload is slow,
    # so we keep the RAW file if it's present AND passes a bin-size
    # sanity check.  Everything else gets regenerated fresh.
    RAW_FILE="${DATA_DIR}/${CELLTYPE}_chr3_raw.${BIN}kb.cool"
    BALANCED_FILE="${DATA_DIR}/${CELLTYPE}_chr3.${BIN}kb.cool"
    MCOOL_FILE="${DATA_DIR}/${CELLTYPE}_chr3.mcool"

    if [ -f "${RAW_FILE}" ]; then
        if ! cooler info "${RAW_FILE}" 2>/dev/null | grep -q "bin-size"; then
            echo "[${CELLTYPE}] Existing RAW file has no bin-size — removing."
            rm -f "${RAW_FILE}"
        fi
    fi
    rm -f "${BALANCED_FILE}" "${MCOOL_FILE}"

    # --- Load into raw cool (matches lldelisle's pipeline) ---
    #
    # Pipeline reference:
    #   https://github.com/lldelisle/Hi-C_reanalysis_Bonev_2017#valid-pair-generation
    #
    # Differences from lldelisle's README command:
    #   • We pass "<chromsizes>:<binsize>" as BINS_PATH instead of a
    #     pre-generated mm10.${BIN}kb.bins file.  Rationale: `cooler
    #     makebins` was renamed to `cooler binnify` in cooler 0.10.x, and
    #     the old name silently fails, leaving a 0-byte bins file that
    #     causes "chromsizes empty" at the cload step.  The shortcut has
    #     no such failure mode.
    #   • --nproc and --max-split are driven by NPROC/MAX_SPLIT env vars
    #     (lldelisle uses -p 1 -s 2).  With NPROC=8 and MAX_SPLIT=8 chr3
    #     is cut into 8 chunks so all 8 workers stay busy during binning —
    #     roughly 4× faster than the -p 4 -s 2 default we had before.
    #
    # Everything else is lldelisle verbatim:
    #   -c2 7 -p2 8  →  chr2, pos2 columns (chr1/pos1 are read from the
    #                   tabix index set by `tabix -s 3 -b 4 -e 4`)
    if [ ! -f "${RAW_FILE}" ]; then
        echo "[${CELLTYPE}] Creating raw cool file at ${BIN} kb (nproc=${NPROC}, max-split=${MAX_SPLIT})..."
        cooler cload tabix \
            --nproc "${NPROC}" \
            -c2 7 -p2 8 -s "${MAX_SPLIT}" \
            --assembly mm10 \
            "${SIZES_FILE}:${BIN}000" "${PAIRS_FILE}" "${RAW_FILE}"
    fi

    # --- Zoomify + balance in one step (all resolutions balanced) ---
    #
    # Idiomatic cooler workflow: go straight from the raw cool to a
    # multi-res mcool, balancing every level as we go.  lldelisle only
    # produces single-resolution .cool because they plot with
    # pyGenomeTracks, but our downstream analysis pipeline expects
    # .mcool.  Resolutions: 1 kb (simulation match), 5 kb (Bonev 2017
    # canonical), 10 kb, 25 kb, 50 kb, 100 kb.
    echo "[${CELLTYPE}] Zoomifying + balancing → multi-res mcool (nproc=${NPROC})..."
    cooler zoomify \
        --nproc "${NPROC}" \
        --balance \
        --balance-args "--cis-only" \
        --resolutions 1000,5000,10000,25000,50000,100000 \
        --out "${MCOOL_FILE}" \
        "${RAW_FILE}"

    # --- Single-res balanced .cool for pyGenomeTracks (optional) ---
    # Extracted from the mcool so it's guaranteed to have bin-size set.
    echo "[${CELLTYPE}] Extracting single-res balanced .cool for plotting..."
    cooler cp "${MCOOL_FILE}::/resolutions/${BIN}000" "${BALANCED_FILE}"

    echo "[${CELLTYPE}] Done: ${MCOOL_FILE}"
done

echo ""
echo "============================================"
echo "Output files in ${DATA_DIR}/"
echo "  single-res .cool : ES_chr3.1kb.cool, CN_chr3.1kb.cool"
echo "  multi-res .mcool : ES_chr3.mcool,   CN_chr3.mcool"
echo ""
echo "Next step — extract Sox2 region:"
echo "  python scripts/process_hic/03_extract_sox2.py \\"
echo "    --mcool-es ${DATA_DIR}/ES_chr3.mcool \\"
echo "    --mcool-cn ${DATA_DIR}/CN_chr3.mcool"
echo "============================================"
