#!/usr/bin/env python
"""Ebf3 Stage 3 figure: regulatory landscape of the simulation window.
Run from the repo root with the ebf3_analysis env:
  python -I scripts/ebf3/plot_regulatory_landscape.py"""
import gzip, os, sys
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
sys.path.insert(0, os.getcwd())
from scripts.analysis import style, annotate_loops as al

CHROM, WS, WE = "chr7", 135_800_000, 137_700_000
BIN = 2000
G = "data/raw/gse96107/GSE96107_{}.bw"
style.apply()
ev = pd.read_csv("results/ebf3/regulatory_evidence.tsv", sep="\t")
ev = ev[ev["in_sim_window"].astype(str) == "True"]
nb = (WE - WS) // BIN
x = WS + BIN * (np.arange(nb) + 0.5)


def track(name):
    v = pyBigWig.open(G.format(name)).stats(CHROM, WS, WE, nBins=nb, type="mean")
    return np.array([np.nan if a is None else a for a in v], float)


def peaks(path, gz=False):
    op = gzip.open if gz else open
    out = []
    with op(path, "rt") as f:
        for l in f:
            r = l.split("\t")
            if r[0] == CHROM and int(r[2]) > WS and int(r[1]) < WE:
                out.append((int(r[1]), int(r[2])))
    return out


genes = []
with gzip.open("data/raw/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
    for l in f:
        r = l.split("\t")
        if len(r) > 8 and r[0] == CHROM and r[2] == "gene" and int(r[4]) > WS and int(r[3]) < WE:
            gt = r[8].split('gene_type "')[1].split('"')[0]
            if gt in ("protein_coding", "lincRNA", "antisense"):
                genes.append((int(r[3]), int(r[4]), r[6], r[8].split('gene_name "')[1].split('"')[0], gt))

rows = [("genes", 0.9)]
marks = [("H3K27ac", "CN"), ("H3K27ac", "NPC"), ("H3K4me1", "CN"), ("H3K4me1", "NPC"),
         ("H3K4me3", "CN"), ("H3K4me3", "NPC")]
rows += [(f"{c} {m}", 0.7) for m, c in marks]
rows += [("IDR peaks", 0.55), ("VISTA / tiers", 0.8), ("E–P", 0.6)]
fig, axes = plt.subplots(len(rows), 1, figsize=(7.2, 7.4), sharex=True,
                         gridspec_kw=dict(height_ratios=[h for _, h in rows], hspace=0.12))
ax = axes[0]
PC_ROW = {"Mgmt": 0, "Ebf3": 1, "9430038I01Rik": 2, "Glrx3": 0}
for k, (s, e, st, n, gt) in enumerate(sorted(genes)):
    pc = gt == "protein_coding"
    y = (PC_ROW.get(n, 0) if pc else k % 3) * 0.3
    ax.add_patch(Rectangle((max(s, WS), y - 0.06), min(e, WE) - max(s, WS), 0.12,
                           color=style.TEXT if pc else style.TEXT2, alpha=1 if pc else 0.45, lw=0))
    if pc:
        ax.text((max(s, WS) + min(e, WE)) / 2, y + 0.1, f"{n} ({st})", ha="center", va="bottom",
                fontsize=6.5, style="italic")
ax.set_ylim(-0.2, 0.95); ax.set_yticks([]); ax.set_ylabel("GENCODE\nvM25", rotation=0, ha="right", va="center")

inputs = {c: track(f"{c}_Input") for c in ("CN", "NPC")}
for ax, (m, c) in zip(axes[1:7], marks):
    r = track(f"{c}_{m}") / inputs[c]
    ax.fill_between(x, 0, np.nan_to_num(r), color=style.CELL[c], lw=0, step="mid")
    ax.axhline(1, color=style.TEXT2, lw=0.5, ls=":")
    ax.set_ylim(0, np.nanpercentile(r, 99.5) * 1.1)
    ax.set_ylabel(f"{c} {m}\n/ input", rotation=0, ha="right", va="center")
    ax.tick_params(axis="y", labelsize=5.5)

ax = axes[7]
pk = [("CN H3K27ac", peaks("data/raw/gse96107/GSE96107_CN_H3K27ac.IDR0.05.filt.narrowPeak.gz", True), style.CELL["CN"]),
      ("NPC H3K27ac", peaks("data/raw/gse96107/GSE96107_NPC_H3K27ac.IDR0.05.filt.narrowPeak.gz", True), style.CELL["NPC"]),
      ("CN tier 3", peaks("data/peaks/CN_enhancers_tier3.narrowPeak"), style.TEXT)]
for k, (lab, p, col) in enumerate(pk):
    for s, e in p:
        ax.add_patch(Rectangle((s, k - 0.3), max(e - s, 3000), 0.6, color=col, lw=0))
    ax.text(WE + 10_000, k, lab, fontsize=6, va="center")
ax.set_ylim(-0.6, 2.6); ax.set_yticks([]); ax.set_ylabel("IDR-0.05\npeaks", rotation=0, ha="right", va="center")

ax = axes[8]
for _, r in ev.iterrows():
    if r["element"] == "Ebf3_promoter":
        continue
    t = str(r["tier"])
    col = style.TEXT if "1" in t else (style.TEXT2 if "2" in t else "#b9b8b3")
    ax.add_patch(Rectangle((r["start"], -0.25), max(r["end"] - r["start"], 4000), 0.5, color=col, lw=0))
    if t != "-" or r["element"] == "hs232":      # label tiered + the borderline hs232 only
        ax.text((r["start"] + r["end"]) / 2, 0.35, f"{r['element']}\n{r['CN_H3K27ac_over_input']:.2f}",
                ha="center", va="bottom", fontsize=5.5, color=style.TEXT)
ax.set_ylim(-0.4, 1.6); ax.set_yticks([])
ax.set_ylabel("VISTA\n(CN K27ac/in)", rotation=0, ha="right", va="center")

ax = axes[9]
al.annotate_arc(ax, al.HS737_BP, al.EBF3_TSS_BP, height=0.8, label="hs737 ↔ Ebf3 promoter, 1.23 Mb")
ax.plot([al.HS737_BP, al.EBF3_TSS_BP], [0, 0], "o", ms=4, color=style.EP_RING)
ax.set_ylim(-0.2, 1.25); ax.set_yticks([])
ax.set_xlim(WS, WE)
ax.set_xticks(np.arange(WS, WE + 1, 200_000)); ax.set_xticklabels([f"{v / 1e6:.1f}" for v in ax.get_xticks()])
ax.set_xlabel("chr7 (Mb, mm10)")
for a in axes:
    a.spines["left"].set_visible(False)
    for xx in (al.HS737_BP, al.EBF3_TSS_BP):
        a.axvline(xx, color=style.EP_RING, lw=0.4, ls="--", zorder=0)
fig.text(0.01, -0.035,
         "ChIP: GSE96107 in vitro NPC / CN bigWigs, mark mean / Input mean per 2 kb bin (descriptive ratio, not a test). "
         "VISTA tier shading: black = tier 1, dark grey = tier 2 (CN H3K27ac/input > 2), light = neither "
         "(unlabelled; values in results/ebf3/regulatory_evidence.tsv).\n"
         "hs737 has NO IDR-0.05 H3K27ac peak in NPC or CN. CN H3K4me1 is globally elevated relative to NPC: "
         "do not compare cell types on H3K4me1. No accessibility (ATAC/DNase) data in GSE96107.",
         fontsize=6, color=style.TEXT2)
style.save(fig, "figures/ebf3/04_regulatory_landscape.png")
print("wrote figures/ebf3/04_regulatory_landscape.png")
