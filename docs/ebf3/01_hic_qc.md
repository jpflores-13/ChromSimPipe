# Ebf3 Stage 1: experimental Hi-C and QC gate

## 1a. Maps built

`scripts/ebf3/download_bonev_hic_chr7.sh` (adapted from
`legacy/scripts/download_bonev_hic.sh`, chr3 → chr7; tabix / `cooler cload
tabix -c2 7 -p2 8` / `cooler zoomify --balance --balance-args "--cis-only"` kept
verbatim). Source: GEO GSE161259 (lldelisle HiCUP 0.6.1 reanalysis of Bonev
2017, mm10). Remote names verified: `GSE161259_ES_chr7…`,
`GSE161259_ncx_NPC_chr7…`, `GSE161259_ncx_CN_chr7…`. Local names keep the
`ncx_` prefix on purpose.

Outputs `data/mcool/{ES,ncx_NPC,ncx_CN}_chr7.mcool`, resolutions 1, 5, 10, 25,
50, 100 kb (1 kb present, as `prepare_hic.py` requires), ICE weights in
`weight`. Wall time 25–28 min per sample (8 CPU). Raw single-resolution
coolers were deleted afterwards (the same counts are in the `.mcool`).

## 1b. The cell-source mismatch (read before interpreting anything)

GSE161259 reanalyses only GSM2533818–21 (ES), GSM2533835–38 (ncx_NPC) and
GSM2533839–42 (ncx_CN, neocortex-sorted). The in vitro NPC (GSM2533822–25)
and in vitro CN (GSM2533826–31) are **not** in it. The CTCF, H3K27ac,
H3K4me3 and RNA used here come from the in vitro NPC / CN. Ebf3 is 30.4 /
78.8 FPKM in in vitro CN but 1.02 / 3.93 in ncx_CN. **Every figure that pairs
ncx_ Hi-C with in vitro ChIP says so in its caption.** The per-GSM `.tar.gz`
archives are a custom `.track` format, not matrices, and were not downloaded.

## 1c. Matched in vitro NPC / CN Hi-C (scoped, not started)

SRA (eutils runinfo, 2026-11-10): SRX2636670–73 (in vitro NPC) = 43 runs,
8.72 billion read pairs, 395 GB `.sra`; SRX2636674–79 (in vitro CN) = 44 runs,
8.71 billion read pairs, 391 GB. All paired-end, HiSeq 2000, 2 × 49 bp.

| Item | Estimate |
|---|---|
| Download | 786 GB of `.sra`, streamed run by run |
| Mapping + pairing (HiCUP/bowtie2 or distiller/bwa) | ~3,000–4,000 core-hours (34.9 billion reads) |
| Wall time | ~1–2 days with 20 parallel 16-core jobs, queue permitting |
| Peak scratch | ~1 TB (FASTQ of ~20 runs in flight) |
| Kept | genome-wide pairs ~100–150 GB per cell type; chr7 `.mcool` ~0.3 GB each |

HiCUP matches the GSE161259 processing (comparability with ES / ncx); distiller
is faster. Rough estimate, not measured. **Awaiting approval.**

## 1d. QC

### Pre-registered rules (written in `scripts/ebf3/hic_qc.py` before any count was seen)

- Endpoint resolution = the finest of 1 / 5 / 10 / 25 kb at which the raw count
  in the single hs737 × Ebf3-promoter pixel is ≥ 30 (Poisson CV ≤ 18 %) in every state.
- E–P O/E = balanced pixel / mean balanced value of the same diagonal **inside
  the simulation window** (chr7:135,800,000–137,700,000). The same definition is
  applied to simulated maps.
- Window edge fails if a 200 kb-window insulation boundary stronger than the
  chr7 75th percentile lies within ±20 kb of it.

hs737 = element midpoint 136,083,776; Ebf3 TSS 137,314,445.

### Library metrics

| | ES | ncx_NPC | ncx_CN |
|---|---|---|---|
| chr7 valid pairs (= 1 kb cooler sum) | 138,707,503 | 147,238,051 | 149,027,070 |
| trans pairs in file | 0 | 0 | 0 |
| cis < 1 kb | 5.4 % | 0.6 % | 0.8 % |
| cis 1–20 kb | 18.0 % | 10.6 % | 13.2 % |
| cis ≥ 20 kb | 76.6 % | 88.8 % | 86.0 % |
| raw contacts in window (10 kb) | 1,360,492 | 1,499,172 | 1,676,558 |
| 1 kb window bins removed by ICE | 6.2 % | 6.8 % | 7.1 % |

**cis/trans: not computable.** The GSE161259 per-chromosome files contain only
chr7–chr7 pairs (0 trans in all three, counted over the full files). ES has ~7×
more < 1 kb pairs than the neural samples, a library difference; P(s) shapes
are compared, not heights.

### hs737 × Ebf3-promoter pixel

| State | 1 kb | 5 kb | 10 kb | 25 kb |
|---|---|---|---|---|
| ES | 0 (O/E 0.00) | 9 (1.66) | **31 (1.35)** | 132 (1.17) |
| ncx_NPC | 0 (0.00) | 24 (2.26) | **92 (2.00)** | 422 (1.80) |
| ncx_CN | 1 (3.93) | 38 (2.71) | **134 (2.24)** | 575 (1.95) |

**Endpoint resolution: 10 kb.** Fixed; not to be changed after simulations
are seen. The 1 kb pixel is empty, so the E–P endpoint is reported at 10 kb even
though the simulation runs at 1 kb. ES clears the rule with 31 counts, the
narrowest margin.

Split-half check (`scripts/ebf3/split_half_coolers.sh`; each pair assigned to
half A or B by a seeded coin flip; `results/ebf3/ep_oe_experimental_halves.json`):

| State | full | half A | half B | |A − B| |
|---|---|---|---|---|
| ES | 1.35 | 1.42 | 1.27 | 0.15 |
| ncx_NPC | 2.00 | 2.06 | 1.96 | 0.10 |
| ncx_CN | 2.24 | 2.10 | 2.38 | 0.28 |

Observation: E–P O/E is higher in both neural states than in ES (1.35 →
2.00 / 2.24), by more than the split-half spread. This is an ensemble contact
frequency difference in neocortex vs ES Hi-C; it is not a mechanism.

### Window edges

| State | Strong boundary near the left edge | Distance from 135,800,000 | Gate (±20 kb) |
|---|---|---|---|
| ES | 135,720,000 (strength 1.20) | 80 kb outside | pass |
| ncx_NPC | 135,760,000 (1.32) | 40 kb outside | pass |
| ncx_CN | 135,770,000–135,780,000 (1.34) | 20 kb outside | pass, **by one bin** |

Right edge 137,700,000: no strong boundary within ±20 kb in any state (nearest:
the 137.47 Mb boundary, 230 kb inside).

**Caveat.** The left edge passes the pre-registered rule but sits on the flank
of a strong boundary (log2 insulation at 135.80 Mb = −0.78 / −0.74 / −0.92 in
ES / ncx_NPC / ncx_CN, near the dip). That boundary and the sites forming it are
just outside the simulated chain, so the model cannot reproduce it; the hs737
pixel is 283 kb from the edge. The brief's fallback window is worse: its left
edge 135,467,000 lies **inside** a strong boundary in ncx_NPC (135,480,000) and
ncx_CN (135,470,000). If the edge is to be moved, 135,700,000 (enclosing the
boundary) is the candidate. Not done: production was already running and the
rule passed.

### Figures

- `figures/ebf3/01_hic_qc.png`: coverage across chr7:135.4–138.1 Mb, P(s) and
  local slope per state, the E–P pixel table.
- `figures/ebf3/02_insulation_window_check.png`: insulation (100 / 200 / 500 kb)
  with both candidate windows marked.
- `figures/ebf3/fig1_experimental_hic.png` / `fig1b_…_zoom.png`: annotated maps
  (see the report).

`check_inputs.sh --configfile config/ChromSimConfig.ebf3.yaml`: 18 passed,
0 warnings, 0 failed (`results/ebf3/check_inputs_ebf3.txt`).
