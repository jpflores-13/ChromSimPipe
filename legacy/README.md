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
| `scripts/extract_ctcf_sites.py` | mm10 CTCF extraction (use `scripts/extract_ctcf_sites.py`) |
| `scripts/plot_hic_quicklook.py` | Sox2 Hi-C quick-look plots |
| `scripts/download_hic.py`, `download_bonev_hic.sh`, `sbatch_bonev_hic.sbatch`, `envs/hic_bonev.yml` | Bonev 2017 (GSE96107) Hi-C download/processing |
| `scripts/run_batch.py` | mESC/neuron batch runner |
| `scripts/build_comparable_merge.py` | Depth-matched merges for the mESC/neuron sweep |
| `scripts/run_analysis_all_fastdyn.py`, `analysis/polymer_dynamics_parallel.py` | An earlier, parallel variant of the polymer-dynamics analysis; `scripts/run_analysis_all.py` replaced it |
| `scripts/migrate_legacy_analysis.py` | Moved analysis files from the pre-2026 per-simulation layout into the flat one |
| `scripts/plot_contact_map_panels.py` | Multi-panel contact-map figures from an older output layout |
| `scripts/check_results.py` | Listed which simulation folders had finished (the workflow's `done/` markers replace it) |
| `scripts/export_conserved_ctcf_bed.py` | Exported CTCF sites shared by every ctcf_type, for the MSD probe selection |
| `analysis/batch_analysis.py` | Single-pass analysis over all result folders, before `scripts/run_analysis_all.py` |

The last six rows are ChromSimPipe's own scripts, moved here in November 2026
because the workflow never calls them; they may not match the current output
layout.

The current pipeline is `ChromSimPipe.sh` → `workflows/ChromSimPipe.snakefile`,
which simulates with the Hansen lab's polysim package. ChromSimPipe's earlier
built-in simulation model (`scripts/lef_dynamics.py`, `scripts/run_simulation.py`,
`analysis/lef_lifetimes.py`) was removed in November 2026; it remains in the git
history.
