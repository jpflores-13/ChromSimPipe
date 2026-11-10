# Ebf3 Stage 0: environment and repository audit

Date: 2026-11-10. Cluster: UNC Longleaf. Repo: ChromSimPipe `main` at
`12d5e80` (origin/main) plus uncommitted local loci edits to
`config/ChromSimConfig.yaml` that predate this project and are untouched.

## Repository state found, and what was done

| Finding | Action |
|---|---|
| Local `main` was 9 commits behind `origin/main` (at `dd621db`). `scripts/sim_features.py`, `scripts/prepare_hic.py`, `check_inputs.sh`, `sweep_1d.sh`, `docs/recipes.qmd`, and the `sticky` / `sticky_cross` / `attraction_energy` / `ctcf_capture: {min, max}` / `loading` / `pausing` condition keys exist only from commit `83bdcd9` onward | Fast-forwarded to `12d5e80` (`git pull --ff-only`) after stashing and restoring the uncommitted loci edit; no conflicts. No SLURM jobs were running from the checkout |
| `envs/cohesin_sim.yml` and `envs/ctcf_extraction.yml` already exist | Verified, not recreated (below) |
| No package can read bigWig / bigBed in either pipeline env | Added `envs/ebf3_analysis.yml` (pyBigWig, UCSC bigBedToBed / bigWigToBedGraph, cooler / cooltools pinned to the sim env's versions). Built at `.conda_envs/ebf3_analysis` on /work |

All files read before writing: README.md, CODE_GUIDE.md,
docs/configuration.qmd, docs/simulation.qmd, docs/recipes.qmd,
docs/running.qmd, configs/parameters.py, scripts/sim_features.py,
scripts/extract_ctcf_sites.py, workflows/ChromSimPipe.snakefile, plus
scripts/run_simulation_shard.py (to confirm how `sticky` / `attraction_energy`
reach polysim).

## Package versions

Simulation env (`.snakemake/conda/0544ae75e3043e68c45ece7022cee694_`, built
from the current `envs/cohesin_sim.yml`, byte-identical):

| Package | Version |
|---|---|
| python | 3.10.22 |
| numpy | 1.26.4 |
| openmm | 8.2 |
| cuda-version | 12.4 |
| cooler | 0.12.1 |
| cooltools | 0.7.1 |
| h5py | 3.16.0 |
| **polysim** | 0.1.0, git `ahansenlab/polysim_bootcamp@7f273ab648e7c606c7212702bb5caa1604810006` (matches the pin) |
| polychrom | 0.1.1, git `open2c/polychrom@11a870cac8a3b168a2e0e11995d2e899c3ceb657` (unpinned in the env file) |

CTCF env (`57ee66192590a6310b4655ded76602c8_`, identical to
`envs/ctcf_extraction.yml`): FIMO (MEME) 5.5.9, bedtools 2.31.1, samtools 1.24.

A second, stale sim env (`da8e6980…`, built from an older env file) has **no
polysim**. `check_inputs.sh` and `sweep_1d.sh` pick the first env whose YAML
mentions `hic2cool` (`head -1`); today that is the correct `0544ae…` by
alphabetical order only. Deleting the stale env would remove the hazard; not
done (not mine to delete).

## GPU node check

`scripts/ebf3/gpu_check.py`, run with the sim env via
`srun -A rc_dphansti_pi -p l40-gpu,volta-gpu,gpu --qos gpu_access --gres=gpu:1`:

```
host g0601.ll.unc.edu
NVIDIA GeForce GTX 1080, 580.126.20
| NVIDIA-SMI 580.126.20   Driver Version: 580.126.20   CUDA Version: 13.0 |
openmm 8.2 platforms: ['Reference', 'CPU', 'CUDA', 'OpenCL']
CUDA context OK, E = 0.011092695407569408 kJ/mol
polysim import OK
```

A CUDA context was created and integrated (kernels compiled with nvrtc) on a
Pascal GPU, which CUDA 13 no longer supports, confirming the 12.4 pin works.
`polysim.extrusion` and `polysim.sim3d` import (pyximport build OK). Shards may
land on l40, volta or gpu (GTX 1080) nodes; wall time will differ by node type.

## SLURM profile

`profiles/longleaf` is correct for this cluster: `sacctmgr` lists account
`rc_dphansti_pi` with QOS `gpu_access, gpu_rtx6k, normal`; partitions
`general`, `l40-gpu`, `volta-gpu`, `gpu` all exist. `simulateShard` requests
1 GPU / 4 CPU / 60 GB / 720 min, as the brief specifies. No new profile.

## Genome

UCSC mm10 `mm10.fa.gz`, md5 `db005b65828db31735f384e4c5787be5` = UCSC
`md5sum.txt`. Primary chromosomes (chr1-19, X, Y, M) extracted to
`data/mm10.primary.fa` and indexed (`.fai` present, so `configs/parameters.py`
now checks chromosome bounds). Sizes identical to UCSC `mm10.chrom.sizes`.
`genome: 'data/mm10.primary.fa'`, `genome_assembly: 'mm10'`. The lab reference
folder `/proj/phanstiel_lab/Reference/mouse` is empty.

## Storage warning

The home filesystem is over its soft quota: 55.2 GB used of 51.2 GB (hard
limit 76.8 GB). Miniforge (and by default its package cache) lives in
`$HOME/miniforge3`. The ebf3 env and its package cache were put on /work
(`CONDA_PKGS_DIRS=/work/users/j/p/jpflores/.conda_pkgs`). Snakemake envs live
in `.snakemake/conda` on /work. Freeing home space before the next env build
is advisable.

## check_inputs.sh

Run after Stage 1a produces the `.mcool` files; result recorded in
`01_hic_qc.md`.

## Config fields / features expected but not found

| Expected (from the brief) | Found |
|---|---|
| A way to give different sticky tiers to different conditions with one `Enhancer_Peaks` column | `sticky.peaks` names a samplesheet column, so each tier needs its own column. Added `Enhancer_Peaks_T2`, `Enhancer_Peaks_T3`; `Enhancer_Peaks` = tier 1 |
| CTCF `signalValue` usable by `ctcf_capture: {min, max}` | `get_ctcf_signals` matches sites to peaks by **name** (col 4). Every GSE96107 peak is named `.`, so all sites would get one signal and `{min: 0.2, max: 0.9}` would silently become a uniform 0.55. Fixed in the data: `data/peaks/{ES,NPC,CN}_CTCF.narrowPeak` give each peak a unique name |
| Orientation for every type the brief wants reported (ES, NPC, CN) | `orientCTCF` only runs for `ctcf_type`s used by a condition (CN). ES and NPC were oriented by calling `scripts/extract_ctcf_sites.py` directly with the rule's exact arguments |
| A summit for a VISTA element | VISTA elements have none; the element midpoint is used (column 10 = midpoint − start) |
| cis/trans ratio from the Stage-1 pairs | GSE161259 per-chromosome files contain only chr7–chr7 pairs (0 trans in all three); not computable from them |
| Accessibility track in GSE96107 | None (no ATAC / DNase); recorded as unavailable in the evidence table |
| ES H3K27ac / H3K4me1 / H3K4me3 peaks or bigWigs in GSE96107 | Not deposited (ES has CTCF, Ring1B, Input only). No ES tier 3 |
