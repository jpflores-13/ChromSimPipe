#! /bin/bash -login
#SBATCH -J ChromSimPipe
#SBATCH -t 10-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p general
#SBATCH --mem=4gb
#SBATCH -o "%x-%j.out"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80

## ChromSimPipe launcher (bagPipes / dietJuicer pattern).
##
## Usage:
##   sbatch ChromSimPipe.sh                         # run the whole workflow
##   sbatch ChromSimPipe.sh --configfile other.yaml # use a different config
##   bash   ChromSimPipe.sh -n                      # dry run: list jobs only
##   bash   ChromSimPipe.sh --conda-create-envs-only
##
## Any extra arguments are passed straight to snakemake.

## Exit if any command fails
set -e

## Config file (overridable with --configfile)
CONFIG="config/ChromSimConfig.yaml"
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2"; shift 2 ;;
        *)            ARGS+=("$1"); shift ;;
    esac
done

## Read a top-level scalar from the config, e.g. cfg pythonVers
cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }

## Load required modules (versions from the config)
module load python/"$(cfg pythonVers)"
module load anaconda/"$(cfg condaVers)"

## Create and activate virtual environment with requirements
[ -d env ] || python3 -m venv env
source env/bin/activate
python3 -m pip install -q -r config/requirements.txt

## Make directory for slurm logs
SLURM_LOGS="$(cfg outputDir)/logs_slurm"
mkdir -p "${SLURM_LOGS}"

## Execute ChromSimPipe snakemake workflow
snakemake -s workflows/ChromSimPipe.snakefile --profile profiles/slurm --configfile "${CONFIG}" --slurm-logdir "${SLURM_LOGS}" -j 100 --max-jobs-per-second 5 --max-status-checks-per-second 0.5 --rerun-incomplete --keep-going -p --latency-wait 500 "${ARGS[@]}"

## Success message
echo "Workflow completed successfully!"
