## Site setup for San Raffaele HPC, sourced by ChromSimPipe.sh and unlock.sh.
## Make `python3` (>= 3.11, for Snakemake 8) and `conda` (>= 24.7.1, which
## Snakemake 8.27 requires) available here.
##
## The cluster has no python or conda modules and its python3 is 3.9, so this uses
## your own Miniforge (one-time install, ~0.5 GB):
##   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
##   bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
## Set CHROMSIM_CONDA_ROOT if you installed it somewhere else. If an administrator
## builds the environments instead, see "Environments built by an administrator"
## in docs/setup.qmd.
source "${CHROMSIM_CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate base
