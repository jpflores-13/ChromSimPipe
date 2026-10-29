# ChromSimPipe — Cohesin Loop-Extrusion Simulation Pipeline

## What this pipeline does

Inside every cell, long DNA molecules are organised into 3D structures by
**cohesin**, a ring-shaped protein that pulls DNA through itself, extruding
a growing loop. **CTCF** proteins act as roadblocks — cohesin stops when it
runs into a CTCF site from the correct direction. The pattern of where
cohesin stops determines which genomic regions contact each other,
measurable by **Hi-C**.

ChromSimPipe is a physics-based polymer simulation framework for testing
mechanistic hypotheses about 3D chromatin organisation. Given CTCF binding
data (CUT&Tag or ChIP-seq peaks) and Hi-C maps as input, it simulates
cohesin loop extrusion under different parameter conditions and compares
the resulting contact maps to your experimental Hi-C.

**Included example.** The default configuration ships with an example from
the STRS project (Flores et al. 2026): HEK293T cells under hyperosmotic
stress (sorbitol treatment), comparing control vs. sorbitol CTCF CUT&Tag
peaks against matched Hi-C maps. Swap in your own data and the pipeline
runs the same way.

**If you want a line-by-line walkthrough of every script, read
[`CODE_GUIDE.md`](CODE_GUIDE.md).**

### How the pieces fit together

```
  ChromSimSamplesheet.txt                      config/ChromSimConfig.yaml
  (one row per sample: metadata,               (locus, conditions, physics,
   CTCF_Peaks, Hi-C)                            genome, motif, envs, outputDir)
          │                                               │
          └───────────────────────┬───────────────────────┘
                                  ▼
                        configs/parameters.py
                (parses both; the only place values are read)
                                  │
                                  ▼
                   workflows/ChromSimPipe.snakefile
                                  │
     convertHic ── getMotif ── orientCTCF ── validateCTCF
                                  │
             simulateShard (condition × replicate × shard, GPU)
                                  │
                     mergeShards (condition × replicate)
                                  │
              analyzeCondition (contact map, P(s), APA,
               insulation, comparison with Hi-C) ── benchmarkSummary
                                  │
                                  ▼
                       output/analysis/{condition}/
```

### Example conditions (STRS use case)

`config/ChromSimConfig.yaml` ships three example conditions comparing two
experimental samples, `condition1` and `condition2` (each a row in
`ChromSimSamplesheet.txt` with its own CTCF peaks and Hi-C map). To define
your own, edit the `conditions:` list there.

| # | Name | Cohesin params | CTCF sites | Compare vs | Question |
|---|------|---------------|------------|------------|----------|
| 1 | `cond1_baseline`             | Gabriele 2022 | Condition 1 | Condition 1 Hi-C | Baseline: can we reproduce condition 1 Hi-C? |
| 2 | `cond2_ctcf-only`            | Gabriele 2022 | Condition 2 | Condition 2 Hi-C | Do condition 2 CTCF sites alone explain condition 2 Hi-C? |
| 3 | `cond2_weak-ctcf`            | Gabriele 2022 (0.5× CTCF capture) | Condition 2 | Condition 2 Hi-C | Weaker CTCF stalling at retained sites? |

### Example locus (STRS use case, hg38, HEK293T)

The default configuration simulates the 1.96 Mb chr1 window used in the STRS
figures. Change the `locus:` block in `config/ChromSimConfig.yaml` to switch
loci (and the `tiling:` block to match: `tile_size` must equal the locus
length in monomers + 2 × `padding`).

| Key | Chrom | Window | Genes |
|-----|-------|--------|-------|
| `chr1_fig1`  | chr1  | 64.21 – 66.17 Mb | JAK1 / AK4 / LEPR |

---

## Prerequisites

- A SLURM cluster with a GPU partition, or one machine with a GPU (see
  [Running on a different cluster](#running-on-a-different-cluster))
- `python3` ≥ 3.11 and `conda` (loaded by your site's `profiles/<site>/env.sh`)
- CTCF peak files (CUT&Tag / ChIP-seq narrowPeak) and `.hic` maps for your samples
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
   experimental sample. `CTCF_Peaks` (narrowPeak/BED) and `Hi-C` (`.hic`)
   are required columns (write `NA` if a sample has no file of that kind);
   every other column is metadata:

    | Project | Cell_Type | Condition  | CTCF_Peaks                  | Hi-C                     |
    |---------|-----------|------------|-----------------------------|--------------------------|
    | STRS    | HEK293T   | condition1 | /path/to/control.narrowPeak | /path/to/control.hic     |
    | STRS    | HEK293T   | condition2 | /path/to/sorbitol.narrowPeak| /path/to/sorbitol.hic    |

3. Edit `config/ChromSimConfig.yaml`. `fileNamesFrom` lists the samplesheet
   columns joined to name each sample; conditions refer to samples by that
   name:

    ```yaml
    samplesheet: 'ChromSimSamplesheet.txt'
    fileNamesFrom: ['Condition']      # → sample names condition1, condition2
    outputDir: 'output'
    genome: '/path/to/genome.fa'      # FIMO motif orientation
    ...
    conditions:
      - name:        cond1_baseline
        ctcf_type:   condition1       # sample whose CTCF peaks to simulate
        compare_hic: condition1       # sample whose Hi-C to compare against
        ...
    ```

4. Pick your compute site with `hpcProfile` in the config (default
   `profiles/longleaf`). For another cluster, see
   [Running on a different cluster](#running-on-a-different-cluster).

5. Submit the workflow:

    ```bash
    sbatch ChromSimPipe.sh
    ```

    `ChromSimPipe.sh` sources the site's `env.sh` (module loads), builds a
    venv from `config/requirements.txt`, and runs Snakemake with the site's
    profile. Each rule becomes its own SLURM job; the conda
    environments in `condaEnvs:` are built by Snakemake on first use.

### Dry run first

```bash
bash ChromSimPipe.sh -n
```

Prints every job and shell command that would be submitted; nothing runs.
Any extra arguments are passed to `snakemake` (e.g.
`bash ChromSimPipe.sh --conda-create-envs-only` builds the environments
without running anything).

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
| `convertHic` | 1 per Hi-C sample (2) | `hic2cool convert` `.hic` → `.cool` at `resolution` |
| `getMotif` | 1 | Downloads `ctcf_motif.jaspar_id` from JASPAR (skipped if `ctcf_motif.meme` is set) |
| `orientCTCF` | 1 per CTCF sample used (2) | FIMO motif strand → oriented CTCF BED for the locus |
| `validateCTCF` | 1 | Checks every BED parses and maps onto the polymer |
| `simulateShard` | conditions × reps × shards (36, GPU) | polychrom 3D loop-extrusion simulation |
| `mergeShards` | conditions × reps (9) | Concatenates shard trajectories |
| `analyzeCondition` | 1 per condition (3) | Contact maps, P(s), insulation, APA, Hi-C comparison |
| `benchmarkSummary` | 1 | Collects per-job run time / memory into one table |

## Running on a different cluster

Everything specific to a compute site lives in one folder under `profiles/`;
the config names it with `hpcProfile`. The science config and samplesheet do
not change when you move.

| Folder | For |
|--------|-----|
| `profiles/longleaf/` | UNC Longleaf (SLURM) |
| `profiles/slurm-template/` | Any other SLURM cluster: copy and fill in the `<...>` values |
| `profiles/local/` | One machine with no scheduler (a workstation, or an interactive GPU node) |

Each folder has two files:

- `env.sh` — sourced by `ChromSimPipe.sh` and `unlock.sh`; makes `python3`
  (≥ 3.11) and `conda` available (`module load ...`, or
  `source ~/miniforge3/etc/profile.d/conda.sh`).
- `config.yaml` — the Snakemake profile: scheduler, account, partitions, and
  each rule's CPUs / memory / GPU.

To set up a new SLURM cluster:

```bash
cp -r profiles/slurm-template profiles/<cluster>
nano profiles/<cluster>/config.yaml   # account, CPU partition, GPU partition (+ QOS if required)
nano profiles/<cluster>/env.sh        # module loads for python3 and conda
# then in config/ChromSimConfig.yaml:  hpcProfile: 'profiles/<cluster>'
sbatch ChromSimPipe.sh
```

`sinfo -s` lists partitions, `sacctmgr show assoc user=$USER` lists your
accounts, and `module avail` shows the available software. Without editing the
config, `--profile profiles/<cluster>` picks a site for a single run. For
`profiles/local`, start with `bash ChromSimPipe.sh` instead of `sbatch`.

---

## Output structure

Everything is written under `outputDir` (default `output/`):

```
output/
├── {runName}_ChromSimSamplesheet.txt   # samplesheet with sample names
├── {runName}_benchmarking.tsv          # wall time / memory per job
├── mcool/{sample}.cool
├── motifs/{jaspar_id}.meme
├── ctcf/ctcf_oriented_{assembly}_{sample}_{chrom}_{start}_{end}.bed
├── sims/
│   ├── {params_name}_ctcf-{sample}_rep{N}_shard{M}/   # per-shard GPU output
│   └── merged_{params_name}_ctcf-{sample}_rep{N}/     # merged per replicate
├── analysis/{condition}/               # contact maps, P(s), APA, figures, metrics
├── logs/{rule}_{wildcards}.out|.err    # each job's stdout / stderr
├── logs_slurm/                         # SLURM job logs
├── benchmarks/{rule}_{wildcards}.tsv
└── done/                               # completion markers Snakemake tracks
```

`{runName}` joins the metadata values shared by every sample (here
`STRS_HEK293T`).

---

## Pipeline files

| File | Purpose |
|------|---------|
| `ChromSimPipe.sh` | sbatch entry point — run this to launch the workflow |
| `unlock.sh` | Release the Snakemake lock after a failed run |
| `ChromSimSamplesheet.txt` | One row per sample: metadata, CTCF peaks, Hi-C |
| `config/ChromSimConfig.yaml` | Locus, conditions, physics, paths, software versions |
| `config/requirements.txt` | Snakemake version for the launcher venv |
| `profiles/<site>/config.yaml` | Scheduler settings and per-rule resources for one compute site |
| `profiles/<site>/env.sh` | Module loads for that site |
| `workflows/ChromSimPipe.snakefile` | The workflow |
| `workflows/utils/benchmarking.py` | Builds `{runName}_benchmarking.tsv` |
| `envs/*.yml` | Conda environments used by the rules |

---

## Cohesin parameters

The default parameter set is derived from **Gabriele et al. 2022** (*Science*)
single-molecule imaging in mESCs. Adjust the `conditions:` list in
`config/ChromSimConfig.yaml` for your biological context.

| Parameter | Default value | Meaning |
|-----------|--------------|---------|
| `lifetime` | 75 steps | ~150 kb processivity at 1 kb/monomer |
| `separation` | 240 monomers | ~8 cohesin rings per 2 Mb locus |
| `ctcf_capture` | 0.125 | 12.5% stalling probability per CTCF encounter |
| `ctcf_release` | 0.0033 | Stalled cohesin lives ~4× longer than free cohesin |

### Tiling trick

The 1,960-monomer locus is simulated as 28 tiled copies on a 70,000-monomer
"chromosome" (2,500 monomers/tile = 1,960 locus + 270 padding each side). This
gives 28× more contact statistics per GPU run. Reference: Yang et al. 2023
*Nat Commun*. All three numbers are set in the `tiling:` block of the config.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Error running conda info` | conda not on `PATH` | Fix the conda line in `profiles/<site>/env.sh` |
| `CUDA not available` | GPU job on wrong partition | Check `simulateShard` under `set-resources:` in `profiles/<site>/config.yaml` |
| `sim_run.engine is 'polychrom' but it cannot be imported` | Sim env missing polychrom | Rebuild the `condaEnvs.sim` environment |
| Snakemake lock error | Previous run crashed | `bash unlock.sh`, then resubmit |
| `Samplesheet ... has no 'CTCF_Peaks' column` | Header typo, or not tab-separated | Headers must be tab-separated and include `CTCF_Peaks` and `Hi-C` |
| `Condition ctcf_type '...' has no CTCF_Peaks` | Condition names a sample that isn't in the samplesheet | Sample names come from `fileNamesFrom`; check the spelling |
| `peaks had no FIMO hit` | `ctcf_motif.unoriented: 'error'` and some peaks lack a motif | Inspect them in `output/logs/orientCTCF_*.err`, or switch to `'drop'` |
| A job failed | — | Its stdout / stderr are in `output/logs/{rule}_{wildcards}.out|.err` |

For the full codebase walkthrough, see [`CODE_GUIDE.md`](CODE_GUIDE.md).

---

## References

### Example dataset (STRS project)
- Flores et al. 2026 — STRS paper. GEO: GSE310051 (Hi-C), GSE310047 (CUT&Tag).

### Loop extrusion model & parameters
- Fudenberg et al., *Cell Reports* **15**:2038–2049 (2016). Canonical loop-extrusion model.
- Banigan et al., *eLife* **9**:e53558 (2020). LEF dynamics and CTCF barriers.
- Gabriele et al., *Science* **376**:496–501 (2022). Live-cell cohesin imaging in mESCs (parameter source).
- Yang et al., *Nat Commun* **14**:1913 (2023). Tiling trick for convergence.

### Hyperosmotic stress chromatin (STRS example context)
- Amat et al., *Genome Res* **29**:18–28 (2019). Hi-C under NaCl osmotic stress.

### Cohesin depletion
- Rao et al., *Cell* **171**:305–320 (2017). Acute cohesin loss abolishes loops.
