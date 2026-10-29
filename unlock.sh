#! /bin/bash -login

## Release the Snakemake lock after a failed or interrupted run, then
## resubmit with sbatch ChromSimPipe.sh.
##
## Usage:
##   bash unlock.sh                          # default config
##   bash unlock.sh --configfile other.yaml

## Exit if any command fails
set -e

CONFIG="config/ChromSimConfig.yaml"
[[ "$1" == "--configfile" ]] && CONFIG="$2"

cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }

## Load required modules
module load python/"$(cfg pythonVers)"

## Activate virtual environment
source env/bin/activate

## Unlock snakemake workflow
snakemake --unlock -s workflows/ChromSimPipe.snakefile --cores 1 --configfile "${CONFIG}"

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun with sbatch ChromSimPipe.sh"
