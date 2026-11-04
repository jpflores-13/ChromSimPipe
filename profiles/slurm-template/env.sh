## Site setup, sourced by ChromSimPipe.sh and unlock.sh.
## Make `python3` (>= 3.11, for Snakemake 8) and `conda` (>= 24.7.1, which
## Snakemake 8.27 requires; check with `conda --version`) available here.
## If the cluster's conda modules are older, use your own Miniforge, as in
## profiles/longleaf/env.sh.
## `module avail python` / `module avail conda` (or miniforge, mamba) show
## what your cluster provides. Examples:
##   module load python/3.12
##   module load miniforge3
## or, with your own install:
##   source ~/miniforge3/etc/profile.d/conda.sh
