# Ebf3 Stage 3: enhancer and promoter sticky sites

Script: `scripts/ebf3/build_regulatory_sites.py` (ebf3_analysis env).
Outputs: `data/peaks/{ES,NPC,CN}_enhancers_tier{1,2}.narrowPeak`,
`data/peaks/{NPC,CN}_enhancers_tier3.narrowPeak`,
`data/peaks/Ebf3_promoter.narrowPeak`, `results/ebf3/regulatory_evidence.tsv`
(+ `_summary.json`). Figure: `figures/ebf3/04_regulatory_landscape.png`.

All files are 10-column narrowPeak on `chr7`. With `span: summit`,
`sim_features.read_peaks` takes the monomer of `start + column 10`.

## Sources and evidence classes

| Class | Source | Status |
|---|---|---|
| annotation | UCSC mm10 `vistaEnhancers.bb` (gbdb, downloaded 2026-11-10); GENCODE vM25 GTF | primary definition of tier 1 / 2 and the promoter |
| peak | GSE96107 IDR-0.05 narrowPeaks: H3K27ac and H3K4me3 for NPC and CN only | ES histone peaks were not deposited |
| signal | GSE96107 bigWigs: mean(mark) / mean(Input) over the element (pyBigWig, exact) | **descriptive ratio, not a statistical test**; bigWig-based, NOT re-quantified from raw reads |
| accessibility | none in GSE96107 (no ATAC / DNase) | **not available**; no proxy used |

`*` VISTA elements have no summit; the element midpoint is used.

## Per-element evidence (window elements, plus hs542 just outside)

| Element | Start | End | Summit* | Monomer | Tier | Evidence class | CN K27ac/in | NPC K27ac/in | CN K4me1/in | NPC K4me1/in | CN K4me3/in | NPC K4me3/in | CN K27ac IDR | NPC K27ac IDR | CN K4me3 IDR | Accessibility |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hs542 | 135,756,833 | 135,757,856 | 135,757,344 | outside | - | annotation(VISTA) | 0.904 | 0.803 | 2.151 | 1.114 | 0.812 | 0.638 | 0 | 0 | 0 | n/a |
| mm1590 | 135,888,195 | 135,890,453 | 135,889,324 | 89 | 2 | annotation(VISTA) | 2.205 | 1.367 | 3.144 | 0.844 | 0.718 | 0.762 | 0 | 0 | 0 | n/a |
| mm2007 | 136,037,914 | 136,040,176 | 136,039,045 | 239 | - | annotation(VISTA) | 1.194 | 0.809 | 1.266 | 0.705 | 1.035 | 1.023 | 0 | 0 | 0 | n/a |
| hs737 | 136,083,204 | 136,084,349 | 136,083,776 | 283 | 1,2 | annotation(VISTA) | 3.146 | 1.238 | 3.735 | 0.866 | 0.554 | 0.546 | 0 | 0 | 0 | n/a |
| mm1574 | 136,389,299 | 136,391,216 | 136,390,257 | 590 | - | annotation(VISTA) | 0.656 | 0.571 | 1.069 | 0.621 | 0.888 | 0.89 | 0 | 0 | 0 | n/a |
| hs699 | 136,490,619 | 136,492,389 | 136,491,504 | 691 | 2 | annotation(VISTA) | 2.446 | 1.056 | 4.452 | 2.088 | 1.034 | 0.999 | 1 | 0 | 0 | n/a |
| hs562 | 136,737,153 | 136,739,340 | 136,738,246 | 938 | - | annotation(VISTA) | 1.461 | 1.035 | 2.63 | 1.075 | 0.935 | 0.846 | 0 | 0 | 0 | n/a |
| hs656 | 137,002,008 | 137,003,274 | 137,002,641 | 1202 | - | annotation(VISTA) | 1.06 | 0.938 | 3.078 | 1.374 | 0.601 | 0.962 | 0 | 0 | 0 | n/a |
| hs696 | 137,005,037 | 137,006,321 | 137,005,679 | 1205 | - | annotation(VISTA) | 1.084 | 0.736 | 2.634 | 0.948 | 0.591 | 0.779 | 0 | 0 | 0 | n/a |
| hs331 | 137,039,206 | 137,041,141 | 137,040,173 | 1240 | - | annotation(VISTA) | 1.616 | 1.327 | 4.029 | 1.586 | 0.731 | 0.649 | 1 | 0 | 0 | n/a |
| hs589 | 137,091,735 | 137,093,087 | 137,092,411 | 1292 | - | annotation(VISTA) | 1.629 | 0.764 | 2.784 | 0.81 | 0.849 | 0.986 | 0 | 0 | 0 | n/a |
| hs2776 | 137,148,030 | 137,149,724 | 137,148,877 | 1348 | - | annotation(VISTA) | 1.104 | 0.564 | 2.954 | 0.568 | 0.551 | 0.584 | 0 | 0 | 0 | n/a |
| hs2775 | 137,241,196 | 137,242,779 | 137,241,987 | 1441 | - | annotation(VISTA) | 1.175 | 0.774 | 3.921 | 1.125 | 0.355 | 0.963 | 0 | 0 | 0 | n/a |
| hs232 | 137,250,059 | 137,251,790 | 137,250,924 | 1450 | - | annotation(VISTA) | 1.794 | 1.427 | 5.516 | 1.124 | 0.667 | 0.811 | 0 | 0 | 0 | n/a |
| hs2593 | 137,258,032 | 137,259,168 | 137,258,600 | 1458 | - | annotation(VISTA) | 0.952 | 0.679 | 2.664 | 0.658 | 0.666 | 0.769 | 0 | 0 | 0 | n/a |
| hs2774 | 137,294,983 | 137,296,432 | 137,295,707 | 1495 | - | annotation(VISTA) | 0.891 | 0.971 | 3.312 | 1.095 | 0.953 | 1.03 | 0 | 0 | 0 | n/a |
| Ebf3_promoter | 137,313,945 | 137,314,945 | 137,314,445 | 1514 | promoter | annotation(GENCODE vM25 TSS) | 7.496 | 0.316 | 1.752 | 0.269 | 10.162 | 0.764 | 1 | 0 | 1 | n/a |

"in" = Input. CN H3K4me1 is globally elevated relative to NPC (median
CN/NPC ratio across the 16 VISTA elements = 2.9), so **do not compare cell types on
H3K4me1**; within-cell-type ranking only.

## Promoter

GENCODE vM25 Ebf3 transcripts start at 137,314,445 / 137,314,414 /
137,313,986 / 137,313,916 (all within 530 bp, in monomers 1513–1514), plus an
internal TSS at **137,267,538** (ENSMUST00000209578, 47 kb into the gene) and a
retained-intron transcript at 137,226,802. NPC and CN H3K4me3 IDR peaks
support only the main cluster (CN summits 137,313,947 / 137,314,146 /
137,314,353; CN H3K27ac summit 137,314,067). None is at 137,267,538. **Used: the
main TSS 137,314,445 only.** Promoter narrowPeak = TSS ± 500 bp, summit = TSS
→ monomer 1514 (`span: peak` would mark 1513–1514). This bin also holds
`anchor_bp: 137314000`.

## Tiers

| Tier | Definition (as written in the brief) | CN members | Notes |
|---|---|---|---|
| 1 | hs737 + Ebf3 promoter | hs737 (monomer 283) | primary test; identical element in every cell type |
| 2 | VISTA elements in the window with CN H3K27ac/input > 2, + promoter | **mm1590 (2.21), hs737 (3.15), hs699 (2.45)** | see discrepancy below |
| 3 | all CN H3K27ac IDR peaks in the window, minus peaks within ±1 kb of any GENCODE TSS of Mgmt, Glrx3, 9430038I01Rik, Ebf3, + promoter | 13 peaks (12 distinct monomers); 8 dropped at those promoters | upper bound, not a prediction |

ES / NPC tier 1 and 2 contain the **same elements** as CN (the set is defined
once from annotation + CN signal); only column 7 differs. NPC tier 3 = NPC's own
H3K27ac IDR peaks with the same exclusion (6 peaks). **No ES tier 3:** GSE96107
has no ES H3K27ac peaks (samplesheet `NA`).

### Discrepancy with the brief's Tier-2 list (needs a decision)

The brief lists tier 2 as hs737 (3.08), hs699 (2.80), hs232 (2.15), "plus mm1590
(2.48) if the window is widened". Re-computed here:

| Element | Brief | exact mean | zoom-level mean | NaN→0 | ±500 bp pad | ±1 kb pad |
|---|---|---|---|---|---|---|
| hs737 | 3.08 | 3.15 | 2.79 | 3.15 | 2.53 | 2.15 |
| hs699 | 2.80 | 2.45 | 2.48 | 2.52 | 1.99 | 1.54 |
| hs232 | 2.15 | **1.79** | 1.88 | 1.84 | 1.64 | 1.87 |
| mm1590 | 2.48 | 2.20 | 2.25 | 2.20 | 1.69 | 1.53 |

1. **mm1590 is inside the 1.9 Mb window** (135,888,195–135,890,453, monomer
   89), not only in the widened one.
2. **hs232 falls below 2 under every averaging tried**, so the stated rule
   excludes it. None of the methods reproduces the brief's numbers.
3. The ratio is fragile: input means over these elements are 0.13–0.19, and
   ±500 bp of padding moves hs699 from 2.45 to 1.99. A hard >2 cutoff on this
   quantity is not robust.

The files currently follow the rule as computed (mm1590, hs737, hs699). The
alternative is the brief's explicit list (hs737, hs699, hs232 [+ mm1590]).
Raw-read re-quantification with a proper background model, or relaxed peak
re-calling, would settle it; neither is done yet.

### Other tier notes

- hs737 has **no IDR-0.05 H3K27ac peak** in NPC or CN (confirmed), so tier 3
  as defined **does not contain hs737**. Tier 3 probes what stickiness among
  acetylated elements does; it does not test the hs737 contact directly.
- CN tier-3 peak 013 (summit 137,309,585, monomer 1509) lies inside Ebf3, 5 kb
  from the promoter bin; it is outside the ±1 kb exclusion and kept.
- CN peaks 004 and 005 overlap and share monomer 427; `read_peaks` collapses
  them to one site (12 distinct monomers).
- Tier-3 peaks near other annotated TSSs (lincRNAs etc.): none within ±1 kb.
  The brief excludes only the four named genes; 23 lincRNA / TEC / antisense /
  pseudogene loci in the window are not excluded.
- E and P summits sharing a 1 kb bin: none, for tiers 1 and 2.
