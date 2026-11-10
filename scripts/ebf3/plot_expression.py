#!/usr/bin/env python
"""
Expression of the protein-coding genes in a window across all five Bonev 2017
RNA-seq states (GSE96107 per-GSM gene-level FPKM tables, 2 replicates each).
Bars = mean, points = replicates. In vitro and neocortex (ncx_) share a colour
per cell type; neocortex bars are hatched.
Usage (repo root, ebf3_analysis env):
  python -I scripts/ebf3/plot_expression.py --region chr7:135800000-137700000 --out figures/ebf3/fig1c_expression.png
"""
import argparse, glob, gzip, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd())
from scripts.analysis import style

STATES = ["ES", "NPC", "CN", "ncx_NPC", "ncx_CN"]
LABEL = {"ES": "ES", "NPC": "NPC (in vitro)", "CN": "CN (in vitro)",
         "ncx_NPC": "NPC (neocortex)", "ncx_CN": "CN (neocortex)"}
COLOR = {"ES": style.CELL["ES"], "NPC": style.CELL["NPC"], "CN": style.CELL["CN"],
         "ncx_NPC": style.CELL["NPC"], "ncx_CN": style.CELL["CN"]}


def load():
    rows = []
    for st in STATES:
        for p in sorted(glob.glob(f"data/raw/gse96107_rna/GSM*_{st}_rep*.txt.gz")):
            rep = p.split("_rep")[-1].split(".")[0]
            d = pd.read_csv(p, sep="\t")
            d["state"], d["rep"] = st, int(rep)
            rows.append(d)
    return pd.concat(rows).rename(columns={"Gene_Symbol": "gene"})


def genes(chrom, ws, we):
    out = []
    with gzip.open("data/raw/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
        for l in f:
            r = l.split("\t")
            if len(r) > 8 and r[0] == chrom and r[2] == "gene" and int(r[4]) > ws and int(r[3]) < we \
                    and 'gene_type "protein_coding"' in r[8]:
                out.append((int(r[3]), r[8].split('gene_name "')[1].split('"')[0]))
    return [g for _, g in sorted(out)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", default="chr7:135800000-137700000")
    ap.add_argument("--out", default="figures/ebf3/fig1c_expression.png")
    ap.add_argument("--highlight", default="Ebf3")
    a = ap.parse_args()
    chrom, r = a.region.split(":")
    ws, we = map(int, r.split("-"))
    style.apply()
    d = load()
    gs = genes(chrom, ws, we)
    d = d[d["gene"].isin(gs)]
    d.to_csv(a.out.replace(".png", ".tsv").replace("figures/", "results/"), sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(1.25 * len(gs) + 1.6, 2.8))
    w = 0.8 / len(STATES)
    for k, st in enumerate(STATES):
        for i, g in enumerate(gs):
            v = d[(d.gene == g) & (d.state == st)]["fpkm"].values
            x = i - 0.4 + w * (k + 0.5)
            ax.bar(x, np.mean(v) + 0.01, w * 0.9, color=COLOR[st], hatch="////" if st.startswith("ncx") else None,
                   edgecolor="white", lw=0.4, label=LABEL[st] if i == 0 else None)
            ax.plot([x] * len(v), v + 0.01, "o", ms=2.2, mfc="white", mec=style.TEXT, mew=0.5)
    ax.set_yscale("log")
    ax.set_ylim(0.01, 1000)
    ax.set_ylabel("FPKM (+0.01)")
    ax.set_xticks(range(len(gs)))
    ax.set_xticklabels(gs, style="italic")
    for t in ax.get_xticklabels():
        if t.get_text() == a.highlight:
            t.set_fontweight("bold")
    ax.axhline(1, color=style.TEXT2, lw=0.5, ls=":")
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.16), fontsize=6.3)
    ax.set_title(f"Expression, protein-coding genes in {a.region} (GSE96107 RNA-seq)", fontsize=8, pad=30)
    fig.text(0.01, -0.08, "Bars: mean FPKM of 2 replicates (points). Hatched: neocortex-sorted (ncx_) samples, the cell "
             "source of the ncx_ Hi-C. Solid NPC / CN: in vitro, the cell source of the CTCF / H3K27ac ChIP. "
             "Gene-level values, not read coverage.", fontsize=5.6, color=style.TEXT2, wrap=True)
    style.save(fig, a.out)
    print("wrote", a.out, "genes:", gs)


if __name__ == "__main__":
    main()
