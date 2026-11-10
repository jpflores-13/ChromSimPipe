#!/usr/bin/env python
"""
Ebf3: ATAC-seq accessibility as an enhancer PROXY.

Source (verified via the ENCODE API, 2026-11-10): ENCSR810HQR, ATAC-seq,
Mus musculus B6NTac;B6NCrl/Lap forebrain tissue, embryo E14.5, mm10.
  ENCFF776MOH  conservative IDR thresholded peaks (reps 1+2), narrowPeak
  ENCFF097JFT  fold change over control bigWig (reps 1+2)
This is bulk embryonic forebrain: a PROXY for the Bonev neocortex (ncx_) cells
(E14.5 cortex) and NOT matched to the in vitro NPC / CN. Label it as a proxy.

Outputs
  data/peaks/ATAC_enhancers.narrowPeak    one row per distinct summit in the window,
      minus summits within +-1 kb of a GENCODE TSS of Mgmt, Glrx3, 9430038I01Rik, Ebf3
      (same exclusion as tier 3); 10 columns, signalValue = ATAC fold enrichment
  results/ebf3/regulatory_evidence.tsv     accessibility columns filled in
Run from the repo root with the ebf3_analysis env.
"""
import gzip, os, sys
import numpy as np
import pandas as pd
import pyBigWig
sys.path.insert(0, os.getcwd())
from scripts.ebf3.build_regulatory_sites import tss_of, EXCLUDE_GENES

CHROM, WS, WE = "chr7", 135_800_000, 137_700_000
PEAKS = "data/raw/encode_atac/ENCFF776MOH.bed.gz"
BW = "data/raw/encode_atac/ENCFF097JFT.bigWig"
FLANK = 1000
LABEL = "ATAC E14.5 forebrain (ENCODE ENCSR810HQR; proxy)"

rows = []
with gzip.open(PEAKS, "rt") as f:
    for l in f:
        r = l.rstrip("\n").split("\t")
        if r[0] != CHROM:
            continue
        s, e, summit = int(r[1]), int(r[2]), int(r[1]) + int(r[9])
        if WS <= summit < WE:
            rows.append((s, e, summit, float(r[6])))
# ENCODE lists one row per summit of a peak region; keep each distinct summit, strongest first
rows = sorted({r[2]: r for r in sorted(rows, key=lambda r: r[3])}.values(), key=lambda r: r[2])
tss = sorted(t for g in EXCLUDE_GENES for t in tss_of(EXCLUDE_GENES, CHROM)[g])
keep = [r for r in rows if not any(r[0] - FLANK <= t < r[1] + FLANK for t in tss)]
with open("data/peaks/ATAC_enhancers.narrowPeak", "w") as f:
    for i, (s, e, summit, fe) in enumerate(keep, 1):
        f.write(f"{CHROM}\t{s}\t{e}\tATAC_{i:03d}\t1000\t.\t{fe:.5f}\t-1\t-1\t{summit - s}\n")
mon = {(r[2] - WS) // 1000 for r in keep}
print(f"ATAC summits in window: {len(rows)}; kept after promoter exclusion: {len(keep)} "
      f"({len(mon)} distinct monomers)")

# accessibility evidence per element
ev = pd.read_csv("results/ebf3/regulatory_evidence.tsv", sep="\t")
bw = pyBigWig.open(BW)
fe_max, peak, peak_fe = [], [], []
for _, r in ev.iterrows():
    v = bw.stats(CHROM, int(r["start"]), int(r["end"]), type="max", exact=True)[0]
    fe_max.append(round(float(v), 2) if v is not None else np.nan)
    hit = [x for x in rows if x[0] < r["end"] and x[1] > r["start"]]
    peak.append(int(bool(hit)))
    peak_fe.append(round(max(x[3] for x in hit), 2) if hit else "")
ev["accessibility_proxy"] = LABEL
ev["ATAC_max_foldchange"] = fe_max
ev["ATAC_IDRpeak"] = peak
ev["ATAC_peak_foldenrich"] = peak_fe
ev.to_csv("results/ebf3/regulatory_evidence.tsv", sep="\t", index=False)
print(ev[["element", "tier", "CN_H3K27ac_over_input", "ATAC_IDRpeak", "ATAC_peak_foldenrich",
          "ATAC_max_foldchange"]].to_string(index=False))
