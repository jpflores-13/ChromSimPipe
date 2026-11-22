#!/usr/bin/env python
"""
Controls-free cohesin-lifetime comparison at the CTCF-occupancy stall setting.

  python scripts/ebf3/plot_lifetime_compare.py <out_stem> <short_dir> <long_dir> [--short-label 75,000] [--long-label 225,000]

<short_dir>, <long_dir>: results/<locus>/panels_occ folders (written by export_panel_matrices.py <locus> occ)
for the short and long cohesin lifetime, on the SAME view. Columns, left to right:
  Hi-C | B (extrusion + CTCF) short | B long | C (sticky E-P, no cohesin; from the long run) | D short | D long
Rows: ES, NPC, CN. Every map: contact / its own mean at 50 kb, one shared log scale (YlGnBu),
the same transform as plot_celltype_rectangles.R. Controls A and A0 are left out.
Writes <out_stem>.png and .pdf.
"""
import argparse, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

BIN = 10_000
RATIO = 1 / 2.3
CELLS = ["ES", "NPC", "CN"]
CELLC = {"ES": "#2a78d6", "NPC": "#eb6834", "CN": "#1baf7a"}
COLC = {"HiC": "#0b0b0b", "B": "#6250d6", "C": "#c48500", "D": "#e34948"}


def dense(path, lo, hi):
    a = np.loadtxt(path, skiprows=1, ndmin=2)
    n = (hi - lo) // BIN
    m = np.full((n, n), np.nan)
    i, j = ((a[:, 0] - lo) // BIN).astype(int), ((a[:, 1] - lo) // BIN).astype(int)
    k = (i >= 0) & (j < n) & (i < n) & (j >= 0)
    m[i[k], j[k]] = a[k, 2]
    m[j[k], i[k]] = a[k, 2]
    return m


def draw(ax, m, lo, ws, we, norm, cmap):
    n = m.shape[0]
    e = np.arange(n + 1) * BIN + lo
    X, Y = np.meshgrid(e, e)                       # corner grid: X = column pos, Y = row pos
    xr, yr = (X + Y) / 2, (X - Y) / 2              # rotate 45 degrees; y = half the separation
    ax.pcolormesh(xr, yr, np.ma.masked_invalid(m), norm=norm, cmap=cmap, rasterized=True, linewidth=0)
    ax.set_xlim(ws, we)
    ax.set_ylim(0, RATIO * (we - ws))
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out"); ap.add_argument("short"); ap.add_argument("long")
    ap.add_argument("--short-label", default="75,000"); ap.add_argument("--long-label", default="225,000")
    A = ap.parse_args()
    ms, ml = (json.load(open(f"{d}/meta.json")) for d in (A.short, A.long))
    assert (ms["start"], ms["end"]) == (ml["start"], ml["end"]), "short and long runs need the same view"
    ws, we = ms["start"], ms["end"]
    marg = int(RATIO * (we - ws)) + BIN
    lo, hi = (ws - marg) // BIN * BIN, -(-(we + marg) // BIN) * BIN
    z0 = min(ms["zrange"][0], ml["zrange"][0]); z1 = max(ms["zrange"][1], ml["zrange"][1])
    norm = LogNorm(vmin=z0, vmax=z1)
    cmap = plt.get_cmap("YlGnBu")

    cols = [("HiC", "Hi-C", A.long, "HiC", None),
            ("B", "B  extrusion + CTCF", A.short, "B_extrusion", A.short_label),
            ("B", "B  extrusion + CTCF", A.long, "B_extrusion", A.long_label),
            ("C", "C  sticky E–P", A.long, "C_sticky_t1", "no cohesin"),
            ("D", "D  extrusion + sticky", A.short, "D_combined_t1", A.short_label),
            ("D", "D  extrusion + sticky", A.long, "D_combined_t1", A.long_label)]
    W = 3.4
    H = W * RATIO
    fig = plt.figure(figsize=(len(cols) * W + 2.3, 3 * H + 1.3))
    gene = ms["gene"]
    tss = gene["end"] if gene["strand"] == "-" else gene["start"]
    marks = [tss] + ([ms["ep"][0]] if ms.get("ep") else [])
    fw, fh = fig.get_size_inches()
    for r, cell in enumerate(CELLS):
        for c, (key, title, d, col, sub) in enumerate(cols):
            ax = fig.add_axes([(0.9 + c * W) / fw, (0.25 + (2 - r) * H) / fh, (W - 0.12) / fw, (H - 0.1) / fh])
            draw(ax, dense(f"{d}/{cell}_{col}.tsv", lo, hi), lo, ws, we, norm, cmap)
            for x in marks:
                ax.axvline(x, color="#c0392b", lw=0.8, alpha=0.7)
            if r == 0:
                ax.set_title(title, fontsize=13, fontweight="bold", color=COLC[key], pad=15)
                if sub:
                    ax.text(0.5, 1.02, ("cohesin lifetime " if sub[0].isdigit() else "") + sub, transform=ax.transAxes,
                            ha="center", va="bottom", fontsize=11, color="#52514e")
            if c == 0:
                fig.text(0.55 / fw, (0.25 + (2 - r) * H + H / 2) / fh, cell, rotation=90, ha="center", va="center",
                         fontsize=15, fontweight="bold", color=CELLC[cell])
    ax = fig.add_axes([(0.9 + len(cols) * W + 0.05) / fw, (0.25 + H) / fh, 0.14 / fw, (2 * H - 0.1) / fh])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=ax)
    cb.set_label("contact / mean at 50 kb (log)", fontsize=10)
    fig.text(0.9 / fw, 0.05 / fh, f"{gene['name']}  {ms['chrom']}:{ws / 1e6:.2f}-{we / 1e6:.2f} Mb   (red lines: "
             f"{'Ebf3 promoter and hs737' if ms.get('ep') else 'Bdnf TSS'})", fontsize=10, color="#52514e")
    for ext in ("png", "pdf"):
        fig.savefig(f"{A.out}.{ext}", dpi=200 if ext == "png" else 150)


if __name__ == "__main__":
    main()
