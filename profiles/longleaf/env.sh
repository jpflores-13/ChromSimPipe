## Site setup for UNC Longleaf, sourced by ChromSimPipe.sh and unlock.sh.
## Make `python3` (>= 3.11, for Snakemake 8) and `conda` available here.
## Load only anaconda: Longleaf's Lmod refuses to load the `python` and
## `anaconda` modules together, and anaconda/2024.02 ships Python 3.11.
module load anaconda/2024.02
