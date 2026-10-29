#! /bin/bash -login

## Release the Snakemake lock after a failed or interrupted run, then
## resubmit with ChromSimPipe.sh.
##
## Usage:
##   bash unlock.sh                                       # default config
##   bash unlock.sh --configfile other.yaml --profile profiles/<site>

## Exit if any command fails
set -e

CONFIG="config/ChromSimConfig.yaml"
PROFILE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        --profile)    PROFILE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

cfg() { sed -nE "s/^$1:[[:space:]]*['\"]?([^'\"#[:space:]]+).*/\1/p" "${CONFIG}" | head -1; }
PROFILE="${PROFILE:-$(cfg hpcProfile)}"

## Load this site's software
source "${PROFILE%/}/env.sh"

## Activate virtual environment
source env/bin/activate

## Unlock snakemake workflow
snakemake --unlock -s workflows/ChromSimPipe.snakefile --cores 1 --configfile "${CONFIG}"

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun with ChromSimPipe.sh"
