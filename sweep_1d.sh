#! /bin/bash -login

## Cheap 1D parameter sweep (cohesin on the CTCF track, no 3D polymer, no GPU).
## Arguments go to scripts/sweep_1d.py; see `bash sweep_1d.sh --help`.
## Needs the conda environments and the oriented CTCF files of a run
## (`sbatch ChromSimPipe.sh --until validateCTCF` makes them). Larger sweeps
## belong on a compute node (srun).
##
## Usage:
##   bash sweep_1d.sh --condition cond1_baseline --lifetime 37500,75000,150000 --out sweep.tsv
##   bash sweep_1d.sh --configfile other.yaml --condition ... --out ...

## Exit if any command fails
set -e

CONFIG="config/ChromSimConfig.yaml"
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIG="$2"; shift 2 ;;
        *)            ARGS+=("$1"); shift ;;
    esac
done

## The simulation environment, from `condaEnvs: sim` in the config
source scripts/sim_env.sh
find_sim_env

CHROMSIM_CONFIG="${CONFIG}" "${SIMENV}/bin/python" scripts/sweep_1d.py "${ARGS[@]}"
