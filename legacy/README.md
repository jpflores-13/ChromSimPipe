# legacy/

Code carried over from the lab's mESC-vs-neuron project
(polymer_sim_tutorials_and_neurons) that ChromSimPipe does not use.
It is kept for reference only and is **not** maintained: paths, loci (Sox2,
mm10 chr3:34–36 Mb), cell types (mESC / neuron) and cluster settings
(`ric.tiget` mambaforge, `cuda` / `workq` partitions, `dgx01`) are hard-coded
for that project and the lab's original cluster.

| Path | What it was for |
|---|---|
| `cluster/` | Pre-Snakemake SLURM submit scripts for the original cluster |
| `notebooks/` | mESC vs neuron comparison notebooks (outputs embed local Windows paths) |
| `scripts/extract_ctcf_sites.py` | mm10 CTCF extraction (use `scripts/extract_ctcf_sites_hg38.py`) |
| `scripts/plot_hic_quicklook.py` | Sox2 Hi-C quick-look plots |
| `scripts/download_hic.py`, `download_bonev_hic.sh`, `sbatch_bonev_hic.sbatch`, `envs/hic_bonev.yml` | Bonev 2017 (GSE96107) Hi-C download/processing |
| `scripts/run_batch.py` | mESC/neuron batch runner |
| `scripts/build_comparable_merge.py` | Depth-matched merges for the mESC/neuron sweep |

The current pipeline is `ChromSimPipe.sh` → `workflows/ChromSimPipe.snakefile`.
