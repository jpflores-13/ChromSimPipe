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

## The simulation environment (it has cooler and hic2cool): Snakemake keeps a
## copy of each environment file next to the environment it built.
SIMENV=$(grep -l hic2cool .snakemake/conda/*.yaml 2>/dev/null | head -1 | sed 's/\.yaml$//')
if [[ -z "${SIMENV}" || ! -x "${SIMENV}/bin/python" ]]; then
    echo "The conda environments are not built yet. On a compute node, run:" >&2
    echo "  bash ChromSimPipe.sh --conda-create-envs-only" >&2
    exit 1
fi

CHROMSIM_CONFIG="${CONFIG}" "${SIMENV}/bin/python" scripts/check_inputs.py
