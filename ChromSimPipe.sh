#! /bin/bash -login
#SBATCH -J ChromSimPipe
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16gb    # building the conda envs on a first run needs > 4 GB
#SBATCH -o "%x-%j.out"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

## ChromSimPipe launcher (bagPipes / dietJuicer pattern).
##
## Usage:
##   sbatch ChromSimPipe.sh                         # SLURM cluster
##   bash   ChromSimPipe.sh                         # profiles/local (no scheduler)
##   bash   ChromSimPipe.sh -n                      # dry run: list jobs only
##   sbatch ChromSimPipe.sh --configfile other.yaml # use a different config
##   sbatch ChromSimPipe.sh --profile profiles/<site>   # override hpcProfile
##
## The compute site is a profile folder: its env.sh loads software, its
## config.yaml tells Snakemake how to run jobs. It is chosen by, in order:
## --profile, $CHROMSIM_PROFILE, then hpcProfile in the config.
## Per-cluster settings that stay out of git (e.g. in ~/.bashrc):
##   export CHROMSIM_PROFILE=profiles/pioneer       # this cluster's profile
##   export CHROMSIM_OUTPUT_ROOT=/scratch/...       # where a relative outputDir goes
## Any other arguments are passed straight to snakemake.

## Exit if any command fails
set -e

## Config file and profile (overridable with --configfile / --profile)
CONFIG="config/ChromSimConfig.yaml"
PROFILE=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        --profile)    PROFILE="$2"; shift 2 ;;
        *)            ARGS+=("$1"); shift ;;
    esac
done

## Read a top-level scalar from the config, e.g. cfg outputDir
cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }
PROFILE="${PROFILE:-${CHROMSIM_PROFILE:-$(cfg hpcProfile)}}"
PROFILE="${PROFILE%/}"

## Output directory, resolved the same way as configs/parameters.py
OUTDIR="$(cfg outputDir)"
OUTDIR="${OUTDIR:-output}"
OUTDIR="${OUTDIR/#\~/$HOME}"
[[ "${OUTDIR}" = /* ]] || OUTDIR="${CHROMSIM_OUTPUT_ROOT:-$PWD}/${OUTDIR}"

## Load this site's software (python3, conda)
source "${PROFILE}/env.sh"

## Create and activate virtual environment with requirements
[ -d env ] || python3 -m venv env
source env/bin/activate
python3 -m pip install -q -r config/requirements.txt

## SLURM job logs (only for profiles that use the SLURM executor)
EXTRA=()
if grep -qE "^executor:[[:space:]]*slurm" "${PROFILE}/config.yaml"; then
    SLURM_LOGS="${OUTDIR}/logs_slurm"
    mkdir -p "${SLURM_LOGS}"
    EXTRA+=(--slurm-logdir "${SLURM_LOGS}" --max-jobs-per-second 5 --max-status-checks-per-second 0.5)
fi

## Execute ChromSimPipe snakemake workflow
snakemake -s workflows/ChromSimPipe.snakefile --profile "${PROFILE}" --configfile "${CONFIG}" --rerun-incomplete --keep-going -p "${EXTRA[@]}" "${ARGS[@]}"

## Success message
echo "Workflow completed successfully!"
