#! /bin/bash -login

## Check the samplesheet's inputs against each other and the config before a
## run: chromosome names, assemblies, Hi-C resolution and normalization,
## peaks in each locus window. Prints a report; nothing is written. Quick
## enough for a login node.
##
## Usage:
##   bash check_inputs.sh                          # default config
##   bash check_inputs.sh --configfile other.yaml
## Needs the conda environments, which the first run (or
## `bash ChromSimPipe.sh --conda-create-envs-only`) builds.

## Exit if any command fails
set -e

CONFIG="config/ChromSimConfig.yaml"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2";  shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

## The simulation environment, from `condaEnvs: sim` in the config
source scripts/sim_env.sh
find_sim_env

CHROMSIM_CONFIG="${CONFIG}" "${SIMENV}/bin/python" scripts/check_inputs.py
