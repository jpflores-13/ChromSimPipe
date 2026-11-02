## Site setup for Case Western Pioneer, sourced by ChromSimPipe.sh and unlock.sh.
## Make `python3` (>= 3.11, for Snakemake 8) and `conda` (>= 24.7.1, which
## Snakemake 8.27 requires) available here.
##
## Pioneer's Miniconda3 modules ship conda 23.x, too old for Snakemake, so this uses
## your own Miniforge instead (one-time install, ~0.5 GB):
##   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
##   bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
## Set CHROMSIM_CONDA_ROOT if you installed it somewhere else.
source "${CHROMSIM_CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate base
