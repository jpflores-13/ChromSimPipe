# ChromSimPipe: Cohesin Loop-Extrusion Simulation Pipeline

## What this pipeline does

Inside every cell, long DNA molecules are organised into 3D structures by
**cohesin**, a ring-shaped protein that pulls DNA through itself, extruding
a growing loop. **CTCF** proteins act as roadblocks: cohesin stops when it
runs into a CTCF site from the correct direction. The pattern of where
cohesin stops determines which genomic regions contact each other,
measurable by **Hi-C**.

ChromSimPipe is a physics-based polymer simulation framework for testing
mechanistic hypotheses about 3D chromatin organisation. Given CTCF binding
data (CUT&Tag or ChIP-seq peaks) and Hi-C maps as input, it simulates
cohesin loop extrusion at one or more loci under different parameter
conditions and compares the resulting contact maps to your experimental
Hi-C. The simulations use the Hansen lab's
[polysim](https://github.com/ahansenlab/polysim_bootcamp) model (polymer +
loop extrusion on polychrom / OpenMM, GPU).

**Included example.** The default configuration ships with an example from
the STRS project (Flores et al. 2026): HEK293T cells under hyperosmotic
stress (sorbitol treatment), comparing control vs. sorbitol CTCF CUT&Tag
peaks against matched Hi-C maps. Swap in your own data and the pipeline
runs the same way.

A script-by-script description of the code is in
[`CODE_GUIDE.md`](CODE_GUIDE.md).

### How the pieces fit together

```
  ChromSimSamplesheet.txt                      config/ChromSimConfig.yaml
  (one row per sample: metadata,               (loci, conditions, polysim,
   CTCF_Peaks, Hi-C, optional                   genome, motif, envs, outputDir)
   CTCF_BigWig / CTCF_Input_BigWig)
          │                                               │
          └───────────────────────┬───────────────────────┘
                                  ▼
                        configs/parameters.py
                (parses both; the only place values are read)
                                  │
                                  ▼
                   workflows/ChromSimPipe.snakefile
                                  │
     convertHic ── getMotif ── preparePolysim      (shared by all loci)
                                  │
     orientCTCF ── validateCTCF ── ctcfOccupancy*        (per locus)
                                  │
     simulateShard (locus × condition × replicate × shard, GPU, polysim)
                                  │
                mergeShards (locus × condition × replicate)
                                  │
         analyzeCondition (locus × condition: contact map, P(s), APA,
               insulation, comparison with Hi-C) ── benchmarkSummary
                                  │
  * only for conditions with ctcf_capture {by: occupancy}
                                  ▼
                 output/loci/{locus}/analysis/{condition}/
```

### Example conditions (STRS use case)

`config/ChromSimConfig.yaml` ships three example conditions, each simulated at
every locus, comparing two experimental samples, `condition1` and `condition2` (each a row in
`ChromSimSamplesheet.txt` with its own CTCF peaks and Hi-C map). To define
your own, edit the `conditions:` list there.

| # | Name | Cohesin params | CTCF sites | Compare vs | Question |
|---|------|---------------|------------|------------|----------|
| 1 | `cond1_baseline`             | polysim bootcamp | Condition 1 | Condition 1 Hi-C | Baseline: can we reproduce condition 1 Hi-C? |
| 2 | `cond2_ctcf-only`            | polysim bootcamp | Condition 2 | Condition 2 Hi-C | Do condition 2 CTCF sites alone explain condition 2 Hi-C? |
| 3 | `cond2_weak-ctcf`            | polysim bootcamp (0.5× CTCF capture) | Condition 2 | Condition 2 Hi-C | Weaker CTCF stalling at retained sites? |

### Loci (STRS use case, hg38, HEK293T)

The default configuration simulates the 1.96 Mb chr1 window used in the STRS
figures. `loci:` in `config/ChromSimConfig.yaml` is a list: add entries
(name, chrom, region_start, region_end) to simulate several windows in one
run; each gets its own folder under `output/loci/`, and the converted Hi-C is
shared. Tiling is worked out per locus from `tiling.padding` and
`tiling.chain_size`.

| Key | Chrom | Window | Genes |
|-----|-------|--------|-------|
| `chr1_fig1`  | chr1  | 64.21 – 66.17 Mb | JAK1 / AK4 / LEPR |

---

## Prerequisites

- A SLURM cluster with a GPU partition, or one machine with a GPU (see
  [Running on a different cluster](#running-on-a-different-cluster))
- `python3` ≥ 3.11 and `conda` ≥ 24.7.1 (loaded by your site's `profiles/<site>/env.sh`;
  on Longleaf and Pioneer, a one-time Miniforge install in `$HOME/miniforge3`:
  see `profiles/longleaf/env.sh`)
- CTCF peak files (CUT&Tag / ChIP-seq narrowPeak) for your samples, and optionally
  their Hi-C maps (`.hic`, `.mcool` or `.cool`) to compare the simulations with
- A genome FASTA with the same chromosome names as the peaks

## Quickstart

The workflow follows the PhanstielLab Snakemake pattern
([bagPipes](https://github.com/PhanstielLab/bagPipes),
[dietJuicer](https://github.com/PhanstielLab/dietJuicer)): edit a samplesheet
and a config, then `sbatch` one launcher.

1. Clone the workflow into a working directory:

    ```bash
    git clone https://github.com/jpflores-13/chromsimpipe.git ChromSimPipe
    cd ChromSimPipe
    ```

2. Edit the tab-separated `ChromSimSamplesheet.txt`. Each row is one
   experimental sample. `CTCF_Peaks` (narrowPeak/BED) and `Hi-C` (`.hic`,
   `.mcool` or `.cool`) are required columns (write `NA` if a sample has no file of that kind).
   `CTCF_BigWig` (CTCF ChIP signal) is optional, needed only by `ctcf_capture: {by: occupancy}`;
   `CTCF_Input_BigWig` (its input) is optional on top of that and, when given, is subtracted as
   background. Every other column is metadata:

    | Project | Cell_Type | Condition  | CTCF_Peaks                | Hi-C                 |
    |---------|-----------|------------|---------------------------|----------------------|
    | STRS    | HEK293T   | condition1 | data/control.narrowPeak   | data/control.hic     |
    | STRS    | HEK293T   | condition2 | data/sorbitol.narrowPeak  | data/sorbitol.hic    |

   Paths can be absolute (`/proj/lab/ctcf/control.narrowPeak`) or relative
   to the repo root. The shipped samplesheet uses `data/<file>` with each
   cluster's copies linked into `data/`, so one samplesheet works on every
   cluster; with a single cluster, full paths are simpler. See
   [Input data](#input-data).

3. Edit `config/ChromSimConfig.yaml` (`loci:` for where to simulate, with
   coordinates as whole numbers in bp on the 1 kb grid; `conditions:` for
   what to compare; see `docs/configuration.qmd` for choosing windows, and
   `docs/recipes.qmd` for ready-made perturbation conditions).
   `fileNamesFrom` lists the samplesheet
   columns joined to name each sample; conditions refer to samples by that
   name:

    ```yaml
    samplesheet: 'ChromSimSamplesheet.txt'
    fileNamesFrom: ['Condition']      # → sample names condition1, condition2
    outputDir: 'output'
    genome: 'data/GRCh38.primary_assembly.genome.fa'   # FIMO motif orientation
    ...
    conditions:
      - name:        cond1_baseline
        ctcf_type:   condition1       # sample whose CTCF peaks to simulate
        compare_hic: condition1       # sample whose Hi-C to compare against
        ...
    ```

4. Pick your compute site. The default is `hpcProfile: 'profiles/longleaf'`
   in the config; on any other cluster set `CHROMSIM_PROFILE` once in
   `~/.bashrc` (see
   [Running on a different cluster](#running-on-a-different-cluster)).

5. Submit the workflow:

    ```bash
    sbatch ChromSimPipe.sh
    ```

    `ChromSimPipe.sh` sources the site's `env.sh` (software setup), builds a
    venv from `config/requirements.txt`, and runs Snakemake with the site's
    profile. Each rule becomes its own SLURM job; the conda
    environments in `condaEnvs:` are built by Snakemake on first use.

### Input data

Linking is optional. The samplesheet and `genome` take absolute paths, so
on one cluster you can point straight at your files and skip this section.

The shipped samplesheet and config name their inputs as `data/<file>`, so
the same files work on several clusters. To use them that way, link each
cluster's copies into `data/` once, from the repo root:

```bash
SRC=/path/to/this/clusters/data
ln -s "$SRC"/STRS_HEK293_eGFP-YAP_CTCF_cont_0h_peaks.narrowPeak    data/
ln -s "$SRC"/STRS_HEK293_eGFP-YAP_CTCF_sorbitol_1h_peaks.narrowPeak data/
ln -s "$SRC"/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic  data/
ln -s "$SRC"/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic data/
ln -s /path/to/GRCh38.primary_assembly.genome.fa                    data/
ls -lL data/
```

Links take no space, and everything in `data/` except its README is
git-ignored. [`data/README.md`](data/README.md) has the details.

**Hi-C formats.** The `Hi-C` column takes `.hic` (converted with
`hic2cool`), `.mcool` (the `resolution` in the config, 1 kb by default, is
copied out of it) or `.cool` (must already be at `resolution`). Each becomes
`<outputDir>/mcool/<sample>.cool`. Balanced `.mcool`/`.cool` files carry a
`weight` column and converted `.hic` files carry `KR`;
`analysis.hic_balance: auto` (the default) picks whichever is there.

### Check the inputs

```bash
bash check_inputs.sh
```

Checks the samplesheet's files against each other and the config without
running anything: chromosome names, assemblies (genome vs. Hi-C chromosome
lengths), the Hi-C resolution and normalization, and CTCF peaks in each
locus window. Each check prints `PASS`, `WARN` or `FAIL`. It needs the conda
environments built (see above).

### Dry run first

```bash
bash ChromSimPipe.sh -n
```

Prints every job and shell command that would be submitted; nothing runs.
Any extra arguments are passed to `snakemake` (e.g.
`bash ChromSimPipe.sh --conda-create-envs-only` builds the environments
without running anything).

### Test run first

`config/ChromSimConfig.test.yaml` is the full config cut down to a cheap
end-to-end check: same samplesheet, locus and physics, but one condition, one
replicate, one shard, and 2,000 saved frames instead of 10,000 per replicate.
It is 10 jobs, one of them on a GPU (about 1-2 h), written to `output_test/` (under `$CHROMSIM_OUTPUT_ROOT` if
set). Its results are a smoke test, not
data.

```bash
bash   ChromSimPipe.sh --configfile config/ChromSimConfig.test.yaml -n                     # dry run
sbatch ChromSimPipe.sh --configfile config/ChromSimConfig.test.yaml --until validateCTCF   # CPU steps only
sbatch ChromSimPipe.sh --configfile config/ChromSimConfig.test.yaml                        # one GPU simulation, end to end
```

Run it on every new cluster before a full run: it exercises environment
building, CTCF orientation, GPU access, and the analysis at a small fraction
of a full run's cost.

### Email notifications

Set an address and the workflow emails you once when it starts, once when it
finishes, and once if any job fails (not once per job):

```bash
export CHROMSIM_EMAIL=you@university.edu   # e.g. in ~/.bashrc
```

The address can also go in `notify: email:` in the config; the environment
variable keeps it out of a shared repo and wins if both are set.
`notify: events:` picks which of `start`, `success`, `error` to send, and
`notify: command:` is the shell command that sends one message (default
`mail -s {subject} {email}`; use e.g. `sendmail` on clusters without `mail`).
Check your cluster can send mail with
`echo test | mail -s test $CHROMSIM_EMAIL`. A failed email never fails the
workflow.

### Using a different config

```bash
sbatch ChromSimPipe.sh --configfile config/my_locus.yaml
```

The file must be a complete config; the Snakefile and every script read it.

### After a failed run

```bash
bash unlock.sh               # release the Snakemake lock
sbatch ChromSimPipe.sh       # --rerun-incomplete is on by default
```

### What runs

| Rule | Jobs (default config) | What it does |
|------|-----------------------|--------------|
| `convertHic` | 1 per Hi-C sample (2) | `.hic` (hic2cool), `.mcool` or `.cool` → one `.cool` at `resolution` |
| `getMotif` | 1 | Downloads `ctcf_motif.jaspar_id` from JASPAR (skipped if `ctcf_motif.meme` is set) |
| `preparePolysim` | 1 | Builds polysim's Cython cohesin model once |
| `orientCTCF` | loci × CTCF samples used (2) | FIMO motif strand → oriented CTCF BED for the locus |
| `validateCTCF` | 1 per locus (1) | Checks every BED parses and maps onto the polymer |
| `ctcfOccupancy` | loci × CTCF samples with `ctcf_capture: {by: occupancy}` (0) | CTCF occupancy per oriented site from the ChIP bigWig, minus input when given (`scripts/ctcf_occupancy.py`) |
| `simulateShard` | loci × conditions × reps × shards (18, GPU) | polysim 3D loop-extrusion simulation |
| `mergeShards` | loci × conditions × reps (9) | Concatenates shard trajectories |
| `analyzeCondition` | loci × conditions (3) | Contact maps, P(s), insulation, APA, Hi-C comparison |
| `benchmarkSummary` | 1 | Collects per-job run time / memory into one table |

## Running on a different cluster

Everything specific to a compute site lives in one folder under `profiles/`.
The config, samplesheet and test config are the same on every cluster; what
differs is set once per cluster, outside git:

| Setting | How | Example (Pioneer) |
|---------|-----|-------------------|
| Compute profile | `export CHROMSIM_PROFILE=...` in `~/.bashrc` | `profiles/pioneer` |
| Where output goes | `export CHROMSIM_OUTPUT_ROOT=...` in `~/.bashrc` (unset: the repo folder) | `/scratch/pioneer/users/$USER/ChromSimPipe` |
| Input files | symbolic links in `data/` ([Input data](#input-data)) | links to `~/ChromSimPipe_data/*` |

The profile is chosen by, in order: `--profile` on the command line,
`$CHROMSIM_PROFILE`, then `hpcProfile` in the config (default
`profiles/longleaf`). A relative `outputDir` (`output`, `output_test`) goes
under `$CHROMSIM_OUTPUT_ROOT`; an absolute one is used as written. Open a new
shell (or `source ~/.bashrc`) after editing `~/.bashrc`, then submit as usual:
the SLURM jobs inherit these variables.

| Folder | For |
|--------|-----|
| `profiles/longleaf/` | UNC Longleaf (SLURM) |
| `profiles/pioneer/` | Case Western Pioneer (SLURM) |
| `profiles/slurm-template/` | Any other SLURM cluster: copy and fill in the `<...>` values |
| `profiles/local/` | One machine with no scheduler (a workstation, or an interactive GPU node) |

Each folder has two files:

- `env.sh`: sourced by `ChromSimPipe.sh` and `unlock.sh`; makes `python3`
  (≥ 3.11) and `conda` available (`module load ...`, or
  `source ~/miniforge3/etc/profile.d/conda.sh`).
- `config.yaml`: the Snakemake profile: scheduler, account, partitions, and
  each rule's CPUs / memory / GPU.

To set up a new SLURM cluster:

```bash
cp -r profiles/slurm-template profiles/<cluster>
nano profiles/<cluster>/config.yaml   # account, CPU partition, GPU partition (+ QOS if required)
nano profiles/<cluster>/env.sh        # module loads for python3 and conda
echo 'export CHROMSIM_PROFILE=profiles/<cluster>' >> ~/.bashrc
source ~/.bashrc
sbatch ChromSimPipe.sh
```

`sinfo -s` lists partitions, `sacctmgr show assoc user=$USER` lists your
accounts, and `module avail` shows the available software.
`--profile profiles/<cluster>` picks a site for a single run. For
`profiles/local`, start with `bash ChromSimPipe.sh` instead of `sbatch`.

---

## Output structure

Everything is written under `outputDir` (default `output/` in the repo folder,
or under `$CHROMSIM_OUTPUT_ROOT` when that is set):

```
output/
├── {runName}_ChromSimSamplesheet.txt   # samplesheet with sample names
├── {runName}_benchmarking.tsv          # wall time / memory per job
├── mcool/{sample}.cool
├── motifs/{jaspar_id}.meme
├── loci/{locus}/
│   ├── ctcf/ctcf_oriented_{assembly}_{sample}.bed
│   ├── sims/
│   │   ├── {params_name}_ctcf-{sample}_rep{N}_shard{M}/   # per-shard GPU output
│   │   └── merged_{params_name}_ctcf-{sample}_rep{N}/     # merged per replicate
│   ├── analysis/{condition}/           # contact maps, P(s), APA, figures, metrics
│   └── done/                           # completion markers Snakemake tracks
├── logs/{locus}/{rule}_{wildcards}.out|.err   # each job's stdout / stderr
├── logs_slurm/                         # SLURM job logs
└── benchmarks/{rule}_{locus}_{wildcards}.tsv
```

`{runName}` joins the metadata values shared by every sample (here
`STRS_HEK293T`). `docs/analysis.qmd` lists the analysis files; `docs/ps.qmd`
explains P(s), its local slope, and how to plot them from `ps_curves.tsv` and
`ps_derivative_curves.tsv`.

---

## Pipeline files

| File | Purpose |
|------|---------|
| `ChromSimPipe.sh` | sbatch entry point; run this to launch the workflow |
| `unlock.sh` | Release the Snakemake lock after a failed run |
| `sweep_1d.sh` | Cheap 1D sweep of cohesin settings on the CTCF track, no GPU (`scripts/sweep_1d.py`) |
| `check_inputs.sh` | Check the samplesheet's files against each other and the config before a run (`scripts/check_inputs.py`) |
| `ChromSimSamplesheet.txt` | One row per sample: metadata, CTCF peaks, Hi-C |
| `config/ChromSimConfig.yaml` | Locus, conditions, physics, paths, notifications, compute site |
| `config/ChromSimConfig.test.yaml` | Small 10-job version of the config for testing a new setup |
| `data/` | Input files for a run, linked in per cluster (git-ignored except `data/README.md`) |
| `config/requirements.txt` | Snakemake version for the launcher venv |
| `profiles/<site>/config.yaml` | Scheduler settings and per-rule resources for one compute site |
| `profiles/<site>/env.sh` | Module loads for that site |
| `workflows/ChromSimPipe.snakefile` | The workflow |
| `workflows/utils/benchmarking.py` | Builds `{runName}_benchmarking.tsv` |
| `envs/*.yml` | Conda environments used by the rules |

---

## Cohesin parameters

The defaults are the Hansen lab polysim bootcamp tutorial values, themselves
based on Gabriele et al. 2022 (*Science*) and Goel et al. 2025. Adjust the `conditions:` list in `config/ChromSimConfig.yaml` for your
biological context, and the shared settings in its `polysim:` block.

| Parameter | Default value | Meaning |
|-----------|--------------|---------|
| `lifetime` | 75,000 LEF timesteps | With `vlef` 0.0025, ~190 kb extruded per leg (~375 kb loops) |
| `separation` | 240 monomers | ~8 cohesin rings per 2 Mb locus |
| `ctcf_capture` | 0.5 | 50% stalling probability per CTCF encounter |
| `polysim: life_boost_stalled` | 4 | Cohesin stalled at CTCF on both legs lives 4× longer |
| `polysim: vlef` | 0.0025 | Probability a leg steps per timestep: sets extrusion speed relative to polymer motion |

Each condition can also use more of polysim's model: no cohesin
(`separation: none`), per-site CTCF strength (`ctcf_capture: {min, max}` by
peak-signal rank, `{by: occupancy}` from the CTCF ChIP bigWig, or
`{by: table, file}` from your own per-site values), other factors as sticky sites (`sticky:`, e.g. YY1, LDB1, condensates),
chromatin compaction (`attraction_energy`), and targeted cohesin loading or
pausing (`loading:`, `pausing:`). `docs/configuration.qmd` lists every option
and `docs/recipes.qmd` has a worked example of each, including calibrating
simulated time to live imaging and a cheap 1D screen (`bash sweep_1d.sh`).

### The simulation model

ChromSimPipe runs the Hansen lab's [polysim](https://github.com/ahansenlab/polysim_bootcamp)
package (Cython cohesin model + bond updater on polychrom / OpenMM), pinned in
`envs/cohesin_sim.yml`. Each oriented CTCF site becomes a stall site with
probability `ctcf_capture` per encounter (+ sites stop leftward-moving legs, −
sites rightward-moving ones). With per-site probabilities from occupancy or a
table, several sites of different peaks on one monomer combine as
1 − ∏(1 − p), the chance that at least one bound CTCF stalls the leg. Cohesin moves one monomer per `1 / vlef` LEF
timesteps, each timestep being `poly_steps_per_block` polymer steps, and every
`saveevery`-th block is saved, after an optional unsaved equilibration
(`initskip`) of at least one cohesin lifetime. ChromSimPipe seeds polysim's C
random number generator (`drand48`) per shard; polysim itself does not, so
unseeded runs would share one random sequence.

On a CTCF-driven control locus (chr7:157.18–158.13 Mb, a loop lost after
sorbitol), polysim matched the control Hi-C pattern (obs/exp Pearson 0.47,
SCC 0.32) with a domain and corner as strong as the Hi-C's. ChromSimPipe's
earlier, simpler model gave the same pattern about 7× too weak, and has been
removed.

### Tiling trick

Each locus is simulated as back-to-back copies on one long chain: the
1,960-monomer chr1 locus becomes 28 copies of 2,500 monomers (1,960 + 270
padding on each side) on a 70,000-monomer chain, giving 28× more contact
statistics per GPU run (Yang et al. 2023 *Nat Commun*). The copy size and
count are worked out per locus from `tiling.padding` and `tiling.chain_size`.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Locus ...: region_start and region_end must be multiples of resolution` / `is past the end of` / `chromosome ... is not in` | A `loci:` window off the 1 kb grid, beyond the chromosome, or with a chromosome name the genome doesn't use | Use the window the message suggests; match the genome's chromosome names (`docs/configuration.qmd`, "Choosing a window") |
| `Error running conda info` | conda not on `PATH` | Fix the conda line in `profiles/<site>/env.sh` |
| `CUDA not available` | GPU job on wrong partition | Check `simulateShard` under `set-resources:` in `profiles/<site>/config.yaml` |
| `polysim cannot be imported`, or `preparePolysim` fails | Sim env built before polysim was added, or no C compiler | Rebuild the `condaEnvs.sim` environment (it now includes polysim and `c-compiler`); read `<outputDir>/logs/preparePolysim.err` |
| `polysim needs frames per shard ... multiple of ...` | `frames_per_replicate / n_shards` × `saveevery` is not a multiple of `blocks_per_updater` | Use the frame count the message suggests |
| Snakemake lock error | Previous run crashed | `bash unlock.sh`, then resubmit |
| `MissingInputException` naming a `data/...` file | That file is not linked into `data/` on this cluster, or the link is broken (or use its absolute path instead) | `ls -lL data/`; add or fix the link ([Input data](#input-data)) |
| Output appears in the repo folder instead of scratch | `CHROMSIM_OUTPUT_ROOT` not set in the shell that ran `sbatch` | `echo $CHROMSIM_OUTPUT_ROOT`; add the export to `~/.bashrc`, `source ~/.bashrc`, resubmit |
| Jobs go to the wrong cluster's partitions or account | `CHROMSIM_PROFILE` not set, so the config default (`profiles/longleaf`) was used | `echo $CHROMSIM_PROFILE`; set it in `~/.bashrc` or pass `--profile` |
| `Samplesheet ... has no 'CTCF_Peaks' column` | Header typo, or not tab-separated | Headers must be tab-separated and include `CTCF_Peaks` and `Hi-C` |
| `Condition ctcf_type '...' has no CTCF_Peaks` | Condition names a sample that isn't in the samplesheet | Sample names come from `fileNamesFrom`; check the spelling |
| `peaks had no FIMO hit` | `ctcf_motif.unoriented: 'error'` and some peaks lack a motif | Inspect them in `output/logs/<locus>/orientCTCF_*.err`, or switch to `'drop'` |
| No emails arrive | `mail` missing or blocked on the cluster | Test `echo test \| mail -s test $CHROMSIM_EMAIL`; change `notify: command:` |
| Launcher job `OUT_OF_MEMORY` on the first run | Conda environments are built inside the launcher job | Build them once on a compute node: `srun ... --mem=16g --pty bash`, then `bash ChromSimPipe.sh --conda-create-envs-only` |
| `nvrtc: error: invalid value for --gpu-architecture` | OpenMM built against CUDA 13, which dropped Pascal (GTX 10xx) and Volta (V100) GPUs | `envs/cohesin_sim.yml` pins `cuda-version=12.4`; rebuild the env (`--conda-create-envs-only` on a compute node) |
| `CUDA_ERROR_UNSUPPORTED_PTX_VERSION` | The env's CUDA is newer than the GPU node's driver supports | Check `nvidia-smi` on a GPU node ("CUDA Version"); `cuda-version` in `envs/cohesin_sim.yml` must not exceed it; rebuild the env |
| `No module named 'h5py'` (or another package) in a rule | An environment built before the env file listed it, or a different `python` came first on `PATH` | Rules now call `"$CONDA_PREFIX/bin/python"`; rebuild the env if the package is genuinely missing |
| `no 'KR' normalization in ...cool` | The `.cool` has a different normalization column | Set `analysis.hic_balance` to `auto` or one it lists (`weight`, `VC`, ...) or `False` for raw counts |
| `analyzeCondition` ends in `TIMEOUT` | Full runs pool 150,000 frames per condition; on a slow node this can exceed the limit | Raise `analyzeCondition: runtime:` in the profile (default 1440 min), `bash unlock.sh`, resubmit; finished pieces are reused (`--skip-existing`) |
| `convertHic` ends in `TIMEOUT` | Large `.hic`: a 3.7 GB map takes ~1.6 h | Raise `convertHic: runtime:` in `profiles/<site>/config.yaml` (default 480 min) |
| `UCSC gene fetch failed ... timed out` | Compute nodes have no internet | Harmless: figures are drawn without the gene track |
| A job failed | - | Its stdout / stderr are in `<outputDir>/logs/{locus}/{rule}_{wildcards}.out|.err` (shared steps: `<outputDir>/logs/`) |

For a description of every script, see [`CODE_GUIDE.md`](CODE_GUIDE.md).

---

## References

### Example dataset (STRS project)
- Flores et al. 2026, STRS paper. GEO: GSE310051 (Hi-C), GSE310047 (CUT&Tag).

### Loop extrusion model & parameters
- Fudenberg et al., *Cell Reports* **15**:2038–2049 (2016). Canonical loop-extrusion model.
- Banigan et al., *eLife* **9**:e53558 (2020). LEF dynamics and CTCF barriers.
- Gabriele et al., *Science* **376**:496–501 (2022). Live-cell cohesin imaging in mESCs (parameter source).
- Yang et al., *Nat Commun* **14**:1913 (2023). Tiling trick for convergence.

### Hyperosmotic stress chromatin (STRS example context)
- Amat et al., *Genome Res* **29**:18–28 (2019). Hi-C under NaCl osmotic stress.

### Cohesin depletion
- Rao et al., *Cell* **171**:305–320 (2017). Acute cohesin loss abolishes loops.
