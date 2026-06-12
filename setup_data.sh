#!/usr/bin/env bash
# ##############################################################################
# filename:    setup_data.sh
# project:     STRSsim
# description: Create conda envs, convert STRS .hic files to .mcool, extract
#              oriented CTCF sites, and validate. Run once from the STRSsim/
#              repo root before launching simulations.
#
# Usage (interactive — recommended for first run):
#   srun --mem=32G --cpus-per-task=4 --pty bash
#   bash setup_data.sh [options]
#
# Usage (batch):
#   sbatch setup_data.sh [options]
#
# Options (all have defaults; override any path without editing the script):
#   --strs-root DIR       Root of the STRS project [/work/users/j/p/jpflores/projects/STRS]
#   --hic-control FILE    Control .hic file (relative to --strs-root/data/processed/hic/maps/)
#   --hic-sorbitol FILE   Sorbitol .hic file
#   --ctcf-control FILE   Control CTCF narrowPeak (relative to --strs-root/.../peaks/)
#   --ctcf-sorbitol FILE  Sorbitol CTCF narrowPeak
#   --genome FILE         Path to hg38 FASTA [/proj/phanstiel_lab/Reference/...]
#   --resolution INT      .mcool resolution in bp [1000]
#   --skip-envs           Skip conda environment creation
# ##############################################################################
#SBATCH --job-name=STRSsim_setup
#SBATCH --output=logs/setup_%j.out
#SBATCH --error=logs/setup_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --partition=general

set -euo pipefail
mkdir -p logs

# =============================================================================
# DEFAULTS — override any of these with the flags above
# =============================================================================
STRS_ROOT="/work/users/j/p/jpflores/projects/STRS"
HIC_CONTROL_NAME="YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
HIC_SORBITOL_NAME="YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"
CTCF_CONTROL_NAME="STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak"
CTCF_SORBITOL_NAME="STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak"
HG38_FASTA="/proj/phanstiel_lab/Reference/human/hg38/fasta/GRCh38.primary_assembly.genome.fa"
RESOLUTION=1000
SKIP_ENVS=0

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --strs-root)      STRS_ROOT="$2";         shift 2 ;;
        --hic-control)    HIC_CONTROL_NAME="$2";  shift 2 ;;
        --hic-sorbitol)   HIC_SORBITOL_NAME="$2"; shift 2 ;;
        --ctcf-control)   CTCF_CONTROL_NAME="$2"; shift 2 ;;
        --ctcf-sorbitol)  CTCF_SORBITOL_NAME="$2"; shift 2 ;;
        --genome)         HG38_FASTA="$2";        shift 2 ;;
        --resolution)     RESOLUTION="$2";        shift 2 ;;
        --skip-envs)      SKIP_ENVS=1;            shift ;;
        -h|--help)        sed -n '2,23p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# Assemble full paths from root + filenames
HIC_CONTROL="${STRS_ROOT}/data/processed/hic/maps/${HIC_CONTROL_NAME}"
HIC_SORBITOL="${STRS_ROOT}/data/processed/hic/maps/${HIC_SORBITOL_NAME}"
CTCF_CONTROL_PEAKS="${STRS_ROOT}/data/processed/cutntag/output/peaks/${CTCF_CONTROL_NAME}"
CTCF_SORBITOL_PEAKS="${STRS_ROOT}/data/processed/cutntag/output/peaks/${CTCF_SORBITOL_NAME}"

# =============================================================================
# STEP -1: Create conda environments (idempotent — safe to re-run)
# =============================================================================
if [[ "${SKIP_ENVS}" != "1" ]]; then
    echo "=== Step -1: Setting up conda environments ==="

    # Load the anaconda module (Longleaf-specific; harmless if already loaded)
    module load anaconda/2024.02 2>/dev/null || true

    # Prefer mamba (ships with anaconda/2024.02; uses ~5-10x less memory than
    # the conda solver and is much faster). Fall back to conda if unavailable.
    # IMPORTANT: use the same runner for both create AND run — mamba stores envs
    # in ~/.local/share/mamba/envs/ which conda run cannot find by name.
    if command -v mamba &>/dev/null; then
        ENV_CREATE="mamba env create"
        ENV_RUN="mamba run"
        echo "  Using mamba solver (faster, lower memory)"
    else
        ENV_CREATE="conda env create"
        ENV_RUN="conda run"
        echo "  Using conda solver (mamba not found)"
    fi

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
    # Still need ENV_RUN so later steps work correctly
    if command -v mamba &>/dev/null; then
        ENV_RUN="mamba run"
    else
        ENV_RUN="conda run"
    fi
fi

# =============================================================================
# READ LOCUS + GENOME ASSEMBLY FROM SimConfig.yaml
# =============================================================================
echo ""
echo "=== Reading locus and settings from config/SimConfig.yaml ==="

# Use python/pyyaml to read the config — avoids a bash YAML parser
_config_vals=$(${ENV_RUN} -n cohesin_sim python - <<'PYEOF'
import yaml, sys
with open("config/SimConfig.yaml") as f:
    cfg = yaml.safe_load(f)
loc = cfg["locus"]
print(f"LOCUS_KEY={loc['name']}")
print(f"LOCUS_CHROM={loc['chrom']}")
print(f"LOCUS_START={loc['region_start']}")
print(f"LOCUS_END={loc['region_end']}")
print(f"GENOME_ASSEMBLY={cfg.get('genome_assembly', 'hg38')}")
PYEOF
)
eval "${_config_vals}"
echo "  Locus: ${LOCUS_KEY}  ${LOCUS_CHROM}:${LOCUS_START}-${LOCUS_END}"
echo "  Assembly: ${GENOME_ASSEMBLY}"

# =============================================================================
# STEP 0: Create directory structure
# =============================================================================
echo "=== Step 0: Creating directory structure ==="
mkdir -p data/hic data/mcool data/genome data/motifs data/ctcf_beds
mkdir -p results/polychrom_3d results/analysis results/figures logs

# =============================================================================
# STEP 1: Symlink .hic files into STRSsim data directory
# =============================================================================
echo ""
echo "=== Step 1: Symlinking .hic files ==="

for hic_var in HIC_CONTROL HIC_SORBITOL; do
    hic_path="${!hic_var}"
    hic_dest="data/hic/$(basename ${hic_path})"
    if [[ ! -f "${hic_path}" ]]; then
        echo "  WARNING: ${hic_path} not found — pass --strs-root or --hic-control/--hic-sorbitol to override"
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
    "data/mcool/control.mcool" \
    "${RESOLUTION}"

convert_hic \
    "data/hic/$(basename ${HIC_SORBITOL})" \
    "data/mcool/sorbitol.mcool" \
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
echo "  1. Verify the locus in config/SimConfig.yaml (locus.name)"
echo "  2. Launch simulations: sbatch SimPipe.sh"
