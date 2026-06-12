#!/usr/bin/env bash
# ##############################################################################
# filename:    setup_data.sh
# project:     ChromSimPipe
# description: Create conda envs, convert .hic files to .mcool, extract
#              oriented CTCF sites, and validate. Run once from the repo
#              root before launching simulations.
#
# Paths are read from config/ChromSimConfig.yaml and ChromSimSamplesheet.txt.
# Edit those files before running this script.
#
# Usage (interactive — recommended for first run):
#   srun --mem=32G --cpus-per-task=4 --pty bash
#   bash setup_data.sh [options]
#
# Usage (batch):
#   sbatch setup_data.sh [options]
#
# Options:
#   --genome FILE         Override genome FASTA path (otherwise read from config)
#   --resolution INT      Override .mcool resolution in bp (otherwise from config)
#   --skip-envs           Skip conda environment creation
# ##############################################################################
#SBATCH --job-name=ChromSimPipe_setup
#SBATCH --output=logs/setup_%j.out
#SBATCH --error=logs/setup_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --partition=general

set -euo pipefail
mkdir -p logs

# =============================================================================
# DEFAULTS — paths are read from config; only these flags override
# =============================================================================
GENOME_OVERRIDE=""
RESOLUTION_OVERRIDE=""
SKIP_ENVS=0

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome)         GENOME_OVERRIDE="$2";      shift 2 ;;
        --resolution)     RESOLUTION_OVERRIDE="$2";  shift 2 ;;
        --skip-envs)      SKIP_ENVS=1;               shift ;;
        -h|--help)        sed -n '2,22p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# =============================================================================
# STEP -1: Create conda environments (idempotent — safe to re-run)
# =============================================================================

# Always load anaconda so conda/mamba are available for all later steps
module load anaconda/2024.02 2>/dev/null || true

# Prefer mamba (ships with anaconda/2024.02; faster, lower memory).
# IMPORTANT: use the same runner for both create AND run — mamba stores envs
# in ~/.local/share/mamba/envs/ which conda run cannot find by name.
if command -v mamba &>/dev/null; then
    ENV_CREATE="mamba env create"
    ENV_RUN="mamba run"
else
    ENV_CREATE="conda env create"
    ENV_RUN="conda run"
fi

if [[ "${SKIP_ENVS}" != "1" ]]; then
    echo "=== Step -1: Setting up conda environments ==="
    echo "  Using: ${ENV_RUN%% *}"

    if ${ENV_CREATE%% *} env list 2>/dev/null | grep -qE "cohesin_sim"; then
        echo "  Already exists: cohesin_sim"
    else
        echo "  Creating: cohesin_sim  (may take 5-10 min)"
        ${ENV_CREATE} -f environment.yml
        echo "  Done: cohesin_sim"
    fi

    if ${ENV_CREATE%% *} env list 2>/dev/null | grep -qE "ctcf_extraction"; then
        echo "  Already exists: ctcf_extraction"
    else
        echo "  Creating: ctcf_extraction"
        ${ENV_CREATE} -f envs/ctcf_extraction.yml
        echo "  Done: ctcf_extraction"
    fi
else
    echo "=== Step -1: Skipping conda env creation (--skip-envs) ==="
fi

# =============================================================================
# READ ALL PATHS FROM ChromSimConfig.yaml (falls back to SimConfig.yaml)
# =============================================================================
echo ""
_CFG_FILE="config/ChromSimConfig.yaml"
[[ ! -f "${_CFG_FILE}" ]] && _CFG_FILE="config/SimConfig.yaml"
echo "=== Reading settings from ${_CFG_FILE} ==="

_config_vals=$(${ENV_RUN} -n cohesin_sim python - "${_CFG_FILE}" <<'PYEOF'
import yaml, sys
cfg_path = sys.argv[1]
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)
loc = cfg["locus"]
samplesheet = cfg.get("samplesheet", "ChromSimSamplesheet.txt")
print(f"LOCUS_KEY={loc['name']}")
print(f"LOCUS_CHROM={loc['chrom']}")
print(f"LOCUS_START={loc['region_start']}")
print(f"LOCUS_END={loc['region_end']}")
print(f"GENOME_ASSEMBLY={cfg.get('genome_assembly', 'hg38')}")
print(f"HIC_CONTROL={cfg.get('hic_control', '')}")
print(f"HIC_SORBITOL={cfg.get('hic_sorbitol', '')}")
print(f"GENOME_FROM_CFG={cfg.get('genome', '')}")
print(f"RESOLUTION={cfg.get('resolution', 1000)}")
print(f"SAMPLESHEET={samplesheet}")
PYEOF
)
eval "${_config_vals}"

# Command-line overrides
[[ -n "${GENOME_OVERRIDE}" ]]     && GENOME_FROM_CFG="${GENOME_OVERRIDE}"
[[ -n "${RESOLUTION_OVERRIDE}" ]] && RESOLUTION="${RESOLUTION_OVERRIDE}"
HG38_FASTA="${GENOME_FROM_CFG}"

echo "  Locus: ${LOCUS_KEY}  ${LOCUS_CHROM}:${LOCUS_START}-${LOCUS_END}"
echo "  Assembly: ${GENOME_ASSEMBLY}"
echo "  Samplesheet: ${SAMPLESHEET}"

# Read CTCF peak paths from samplesheet
CTCF_CONTROL_PEAKS=$(awk -F'\t' 'NR>1 && $1=="control"  {print $2; exit}' "${SAMPLESHEET}" 2>/dev/null || echo "")
CTCF_SORBITOL_PEAKS=$(awk -F'\t' 'NR>1 && $1=="sorbitol" {print $2; exit}' "${SAMPLESHEET}" 2>/dev/null || echo "")

# =============================================================================
# STEP 0: Create directory structure
# =============================================================================
echo "=== Step 0: Creating directory structure ==="
mkdir -p data/hic data/mcool data/genome data/motifs data/ctcf_beds
mkdir -p results/polychrom_3d results/analysis results/figures logs

# =============================================================================
# STEP 1: Symlink .hic files into data directory
# =============================================================================
echo ""
echo "=== Step 1: Symlinking .hic files ==="

for hic_var in HIC_CONTROL HIC_SORBITOL; do
    hic_path="${!hic_var}"
    hic_dest="data/hic/$(basename ${hic_path})"
    if [[ ! -f "${hic_path}" ]]; then
        echo "  WARNING: ${hic_path} not found — edit hic_control/hic_sorbitol in ${_CFG_FILE}"
    elif [[ ! -e "${hic_dest}" ]]; then
        ln -s "${hic_path}" "${hic_dest}"
        echo "  Linked: ${hic_dest}"
    else
        echo "  Already exists: ${hic_dest}"
    fi
done

# =============================================================================
# STEP 2: Convert .hic → .mcool
# =============================================================================
echo ""
echo "=== Step 2: Converting .hic → .mcool (resolution: ${RESOLUTION} bp) ==="
echo "    This takes ~5-15 min per file on Longleaf. Run interactively or"
echo "    submit as a CPU job: srun --mem=32G --cpus-per-task=4 --pty bash"

convert_hic() {
    local hic_in="$1"
    local mcool_out="$2"
    local res="$3"

    if [[ ! -f "${hic_in}" ]]; then
        echo "  SKIP: ${hic_in} not found"
        return
    fi

    if [[ -f "${mcool_out}" ]]; then
        echo "  Already exists: ${mcool_out}"
        return
    fi

    echo "  Converting: $(basename ${hic_in}) → $(basename ${mcool_out})"
    ${ENV_RUN} -n cohesin_sim hic2cool convert "${hic_in}" "${mcool_out}" -r "${res}"
    echo "  Done: ${mcool_out}"
}

convert_hic \
    "data/hic/$(basename ${HIC_CONTROL})" \
    "data/mcool/control.cool" \
    "${RESOLUTION}"

convert_hic \
    "data/hic/$(basename ${HIC_SORBITOL})" \
    "data/mcool/sorbitol.cool" \
    "${RESOLUTION}"

# =============================================================================
# STEP 3: Symlink CTCF peak files
# =============================================================================
echo ""
echo "=== Step 3: Symlinking CTCF CUT&Tag peak files ==="

for peaks_var in CTCF_CONTROL_PEAKS CTCF_SORBITOL_PEAKS; do
    peaks_path="${!peaks_var}"
    peaks_dest="data/ctcf_beds/$(basename ${peaks_path})"
    if [[ ! -f "${peaks_path}" ]]; then
        echo "  WARNING: ${peaks_path} not found"
    elif [[ ! -e "${peaks_dest}" ]]; then
        ln -s "${peaks_path}" "${peaks_dest}"
        echo "  Linked: ${peaks_dest}"
    else
        echo "  Already exists: ${peaks_dest}"
    fi
done

# =============================================================================
# STEP 4: Extract oriented CTCF sites for each locus
# =============================================================================
echo ""
echo "=== Step 4: Extracting oriented CTCF sites ==="
if [[ -f "${HG38_FASTA}" ]]; then
    echo "    Using genome FASTA: ${HG38_FASTA}"
else
    echo "    WARNING: hg38 FASTA not found at ${HG38_FASTA}"
    echo "    The extraction script will download ~900 MB — expect 20-30 min extra."
fi

ctcf_control_peaks="data/ctcf_beds/$(basename ${CTCF_CONTROL_PEAKS})"
ctcf_sorbitol_peaks="data/ctcf_beds/$(basename ${CTCF_SORBITOL_PEAKS})"

chrom="${LOCUS_CHROM}"
start="${LOCUS_START}"
end="${LOCUS_END}"
region="${chrom}:${start}-${end}"

for cond in control sorbitol; do
    if [[ "${cond}" == "control" ]]; then
        peaks_file="${ctcf_control_peaks}"
    else
        peaks_file="${ctcf_sorbitol_peaks}"
    fi

    out_bed="data/ctcf_beds/ctcf_oriented_${GENOME_ASSEMBLY}_${cond}_${chrom}_${start}_${end}.bed"

    if [[ -f "${out_bed}" ]]; then
        echo "  Already exists: ${out_bed}"
        continue
    fi

    if [[ ! -f "${peaks_file}" ]]; then
        echo "  SKIP: peaks file not found: ${peaks_file}"
        continue
    fi

    echo "  Extracting: ${LOCUS_KEY} / ${cond} → $(basename ${out_bed})"

    genome_args=()
    [[ -f "${HG38_FASTA}" ]] && genome_args=(--genome "${HG38_FASTA}")

    ${ENV_RUN} -n ctcf_extraction python scripts/extract_ctcf_sites_hg38.py \
        --source bed \
        --bed "${peaks_file}" \
        --region "${region}" \
        "${genome_args[@]+"${genome_args[@]}"}" \
        --output "${out_bed}"

    echo "  Done: ${out_bed}"
done

# =============================================================================
# STEP 5 (removed): BED files now live directly in data/ctcf_beds/ —
# parameters.py reads them there, no symlink needed.
# =============================================================================

# =============================================================================
# STEP 6: Validate CTCF BED files
# =============================================================================
echo ""
echo "=== Step 6: Validating oriented CTCF BED files ==="
if ${ENV_RUN} -n cohesin_sim python scripts/validate_ctcf.py; then
    echo "  Validation passed."
else
    echo "  WARNING: Validation reported issues — check output above before running simulations."
fi

echo ""
echo "=== Setup complete ==="
echo ""
echo "Next steps:"
echo "  1. Verify the locus in config/ChromSimConfig.yaml (locus.name)"
echo "  2. Launch simulations: sbatch ChromSimPipe.sh"
