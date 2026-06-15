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
[`CODE_GUIDE.md`](CODE_GUIDE.md) or [`CODE_GUIDE.pdf`](CODE_GUIDE.pdf).**

### How the pieces fit together

```
  CTCF peaks (CUT&Tag / ChIP-seq)      Experimental Hi-C (.hic maps)
    condition A + condition B              condition A + condition B
          │                                       │
          ▼                                       ▼
  data/ctcf_beds/*.bed             data/mcool/*.cool
   (oriented CTCF sites,            (1 kb contact maps; converted by hic2cool)
    from FIMO JASPAR MA0139.1)
          │                                       │
          └──────────────┬────────────────────────┘
                         ▼
                 configs/parameters.py
           (N conditions × active locus)
                         │
                         ▼
         GPU polymer simulation (polychrom / OpenMM)
         SLURM jobs: N conditions × reps × shards
                         │
                         ▼
              results/polychrom_3d/merged_*/
                         │
                         ▼
              scripts/run_analysis_all.py
           (contact map, P(s), APA, MSD,
            experimental Hi-C comparison)
                         │
                         ▼
                results/analysis/*.npy, *.png
```

### Example conditions (STRS use case)

The default `configs/parameters.py` defines five conditions for the STRS
hyperosmotic-stress example (HEK293T cells, hg38). To define your own
conditions, edit `SIMULATION_CONDITIONS` in `configs/parameters.py`.

| # | Name | Cohesin params | CTCF sites | Compare vs | Question |
|---|------|---------------|------------|------------|----------|
| 1 | `control_ctcf-control`          | Gabriele 2022 | Control CUT&Tag  | Control Hi-C  | Baseline: can we reproduce untreated Hi-C? |
| 2 | `control_ctcf-sorbitol`         | Gabriele 2022 | Sorbitol CUT&Tag | Sorbitol Hi-C | Does CTCF loss alone explain sorbitol Hi-C? |
| 3 | `weak_ctcf-sorbitol`            | Gabriele 2022 (0.5× CTCF capture) | Sorbitol CUT&Tag | Sorbitol Hi-C | Weaker CTCF stalling at retained sites? |
| 4 | `sorbitol_promoter-stall`       | Gabriele 2022 | Sorbitol + promoter (SP/KLF) | Sorbitol Hi-C | Key test: promoter-anchored cohesin? |
| 5 | `sorbitol_promoter-stall_long`  | Gabriele 2022 (2× processivity) | Sorbitol + promoter | Sorbitol Hi-C | Same + longer processivity? |

### Example loci (STRS use case, hg38, HEK293T)

The default configuration ships four ~2 Mb loci from the STRS paper figures.
Change `ACTIVE_LOCUS` in `configs/parameters.py` to switch loci, or add your own.

| Key | Chrom | Window | Genes |
|-----|-------|--------|-------|
| `chr1_fig1`  | chr1  | 64.2 – 66.2 Mb | JAK1 / AK4 / LEPR |
| `chr4_fig1`  | chr4  | 1.7 – 3.7 Mb   | RNF4 / ADD1 / HTT |
| `chr6_fig1`  | chr6  | 124.1 – 126.1 Mb | NKAIN2 / RNF217 |
| `chr16_sox8` | chr16 | 0.4 – 2.3 Mb   | SOX8 |

---

## Prerequisites

- Longleaf HPC account with access to `rc_dphansti_pi` allocation
- `.hic` Hi-C maps and CTCF peak `.narrowPeak` files for your experiment
  (the defaults point to `~/projects/STRS` for the included example)
- `module load anaconda/2024.02` available (already on Longleaf)
- GPU partition access (`general` + `gpu`)

---

## Setup (run once)

### Step 0 — Clone and configure

```bash
cd ~/projects/ChromSimPipe   # or wherever you cloned it

# 1. Edit config/ChromSimConfig.yaml to point at your Hi-C and CTCF data.
nano config/ChromSimConfig.yaml

# 2. Edit ChromSimSamplesheet.txt with your CTCF peak file paths.
nano ChromSimSamplesheet.txt

# 3. Edit profiles/slurm/config.yaml to set your SLURM account.
nano profiles/slurm/config.yaml   # change slurm_account

# 4. Optionally change the locus in ChromSimConfig.yaml (default: chr1_fig1)
```

### Step 1 — Build conda environments and prepare data

```bash
bash setup_data.sh
```

This script reads all paths from `config/ChromSimConfig.yaml` and
`ChromSimSamplesheet.txt`, then:
1. Creates `cohesin_sim` conda env (polychrom, cooler, hic2cool, numpy, scipy, matplotlib)
2. Creates `ctcf_extraction` conda env (MEME/FIMO, bedtools, samtools)
3. Converts input `.hic` files → `.cool` at 1 kb resolution (single-resolution cooler)
4. Extracts and orients CTCF sites from CUT&Tag peaks via FIMO
5. Validates the resulting BED files

Takes ~20–30 minutes total (mostly conda + FIMO).

---

## Running the pipeline

```bash
sbatch ChromSimPipe.sh
```

That's it. `ChromSimPipe.sh` submits itself as a SLURM job, sets up a Python venv
with Snakemake, and launches the full DAG via `snakemake-executor-plugin-slurm`.
Each rule becomes its own SLURM job with the right resources (GPU for
simulation, CPU for everything else).

### What happens

```
convert_hic (×2)
    → extract_ctcf (×2)
        → validate_ctcf
            → simulate_shard (×60: 5 cond × 3 rep × 4 shards, GPU)
                → merge_shards (×15: 5 cond × 3 rep)
                    → analyze_all (×5: one per condition)
```

### Dry run first

```bash
bash ChromSimPipe.sh --dry-run
```

Prints all rules and shell commands that would be submitted — nothing runs.

### Changing locus

```bash
# Edit configs/parameters.py:
ACTIVE_LOCUS = "chr4_fig1"   # or chr6_fig1, chr16_sox8

# Re-run setup for the new locus CTCF beds, then relaunch:
bash setup_data.sh --skip-envs   # skips conda env creation
sbatch ChromSimPipe.sh
```

### After a failed run

```bash
bash unlock.sh   # release Snakemake lock
sbatch ChromSimPipe.sh --rerun-incomplete   # (--rerun-incomplete is the default)
```

---

## Output structure

```
data/
├── hic/                          # symlinks to input .hic files
├── mcool/
│   ├── condition_A.cool          # converted by convert_hic rule (single-resolution)
│   └── condition_B.cool
└── ctcf_beds/
    ├── ctcf_oriented_hg38_condA_{chrom}_{start}_{end}.bed
    └── ctcf_oriented_hg38_condB_{chrom}_{start}_{end}.bed

results/
├── polychrom_3d/
│   ├── {params}_ctcf-{type}_rep{N}_shard{M}/   # raw GPU output (per shard)
│   │   ├── params.json
│   │   └── blocks_*.h5
│   └── merged_{params}_ctcf-{type}_rep{N}/      # merged by merge_shards rule
│       ├── params.json
│       └── blocks_*.h5
│
├── analysis/{condition}/         # per-condition analysis outputs
│   ├── *_contact_map.npy
│   ├── *_ps_curve.npz
│   ├── *_sim_vs_exp_map.png
│   ├── *_msd_*.png / .npz
│   └── ...
│
├── figures/                      # contact map panels
│
└── .sentinels/                   # Snakemake tracking files (safe to ignore)
```

---

## Snakemake pipeline files

| File | Purpose |
|------|---------|
| `ChromSimPipe.sh` | sbatch entry point — run this to launch the pipeline |
| `unlock.sh` | Release Snakemake lock after a failed run |
| `workflows/ChromSimPipe.snakefile` | Full DAG definition (6 rules) |
| `profiles/slurm/config.yaml` | SLURM executor defaults (account, partitions, memory) |
| `config/ChromSimConfig.yaml` | Path + condition configuration (edit before first run) |
| `ChromSimSamplesheet.txt` | CTCF peak file paths (one row per condition) |

---

## Cohesin parameters

The default parameter set is derived from **Gabriele et al. 2022** (*Science*)
single-molecule imaging in mESCs. Adjust `SIMULATION_CONDITIONS` in
`configs/parameters.py` for your biological context.

| Parameter | Default value | Meaning |
|-----------|--------------|---------|
| `lifetime` | 75 steps | ~150 kb processivity at 1 kb/monomer |
| `separation` | 240 monomers | ~8 cohesin rings per 2 Mb locus |
| `ctcf_capture` | 0.125 | 12.5% stalling probability per CTCF encounter |
| `ctcf_release` | 0.0033 | Stalled cohesin lives ~4× longer than free cohesin |

### Tiling trick

The 2 Mb locus (~2000 monomers) is simulated as 28 tiled copies on a 70,000-monomer
"chromosome" (2500 monomers/tile = 2000 locus + 250 padding each side). This gives
28× more contact statistics per GPU run. Reference: Yang et al. 2023 *Nat Commun*.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `ImportError: polychrom` | `cohesin_sim` env not built | Run `bash setup_data.sh` |
| `CUDA not available` | Wrong SLURM partition | GPU jobs need `--partition=gpu --gpus=1` (set automatically by Snakemake) |
| Snakemake lock error | Previous run crashed | Run `bash unlock.sh` then resubmit |
| CTCF BED file not found | `setup_data.sh` not run or CTCF extraction failed | Check `logs/setup_data*.log` |
| No peaks file in samplesheet | Wrong data path | Edit `ChromSimSamplesheet.txt` with correct absolute paths |
| `KeyError: 'CTCF_Type'` in samplesheet | TSV format issue | File must be tab-separated with headers `CTCF_Type` and `CTCF_Peaks_Path` |
| `hic2cool` fails | `.hic` file not found | Check `config/ChromSimConfig.yaml` paths for `hic_control`/`hic_sorbitol` |
| Analysis fails with "No results found" | Merge step incomplete | Check `results/polychrom_3d/` for `merged_*` directories |

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
