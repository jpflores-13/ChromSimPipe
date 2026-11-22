#!/usr/bin/env python
"""
Ebf3 / hs737: experimental Hi-C beside the five models, one row per cell type.

Rows ES, NPC, CN. Columns: Hi-C, then A, A0, B, C, D simulated with that cell
type's CTCF barriers (and, for C/D, its tier-1 enhancer + Ebf3 promoter sites).
View: hs737 - 200 kb to Ebf3 TSS + 200 kb (whole Ebf3 gene inside), 10 kb bins.

Sim maps: output_ebf3/loci/ebf3_hs737/analysis/<condition>/<params>_ctcf-<CT>_30000blk_pooled_contact_map.npy
(1 kb over the simulated window chr7:135.8-137.7 Mb, 3 replicates pooled).
CN conditions carry no suffix (A_polymer ...); ES / NPC ones end in _ES / _NPC.
A missing map is drawn as a "pending" panel, so the figure can be rebuilt as runs finish.

Each map is divided by its own mean contact at 50 kb (sim and Hi-C units differ) and
drawn on one shared log scale. E-P O/E is the pre-registered endpoint (10 kb hs737 x
Ebf3-promoter pixel / mean of its diagonal over the full simulated window).

  figures/ebf3/fig5_celltype_model_panels.{png,pdf}
Run from the repo root with the ebf3_analysis env.
"""
import os, sys
import numpy as np
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sys.path.insert(0, os.getcwd())
from scripts.analysis import style, annotate_loops as al

CHROM = "chr7"
SIM_S, SIM_E = 135_800_000, 137_700_000        # simulated window (config locus)
RES, BIN = 1000, 10_000
F = BIN // RES
EBF3 = (137_193_673, 137_314_445)               # GENCODE vM25 gene body, minus strand
BUFFER = 200_000
WS = (al.HS737_BP - BUFFER) // BIN * BIN        # 135,880,000
WE = -(-(al.EBF3_TSS_BP + BUFFER) // BIN) * BIN  # 137,520,000
GENES = [(136_894_614, 137_130_266, "+", "Mgmt"), (EBF3[0], EBF3[1], "-", "Ebf3"),
         (137_375_574, 137_410_762, "-", ""), (137_437_614, 137_468_594, "+", "Glrx3")]
ANALYSIS = "output_ebf3/loci/ebf3_hs737/analysis"
CELLS = [("ES", "ES", "data/mcool/ES_chr7.mcool"),
         ("NPC", "ncx_NPC", "data/mcool/ncx_NPC_chr7.mcool"),
         ("CN", "ncx_CN", "data/mcool/ncx_CN_chr7.mcool")]
MODELS = [("A_polymer", "polymer"), ("A0_polymer_attr", "polymer_attr"), ("B_extrusion", "extrusion"),
          ("C_sticky_t1", "sticky_t1"), ("D_combined_t1", "combined_t1")]
OUT = "figures/ebf3/fig5_celltype_model_panels.png"
style.apply()


def coarsen(m):
    n = m.shape[0] // F
    m = np.nan_to_num(m[: n * F, : n * F])
    return m.reshape(n, F, n, F).sum(axis=(1, 3))


def ep_oe(m10):
    i, j = (al.HS737_BP - SIM_S) // BIN, (al.EBF3_TSS_BP - SIM_S) // BIN
    return float(m10[i, j] / np.nanmean(np.diagonal(m10, j - i)))


def crop(m10):
    a, b = (WS - SIM_S) // BIN, (WE - SIM_S) // BIN
    return m10[a:b, a:b]


def norm50(m):
    return m / np.nanmean(np.diagonal(m, 50_000 // BIN))


def sim_map(cell, cond, params):
    name = cond if cell == "CN" else f"{cond}_{cell}"
    p = f"{ANALYSIS}/{name}/{params}_ctcf-{cell}_30000blk_pooled_contact_map.npy"
    return coarsen(np.load(p)) if os.path.exists(p) else None


def main():
    grid = {}
    for cell, hic_label, mcool in CELLS:
        hic = cooler.Cooler(f"{mcool}::resolutions/{BIN}").matrix(balance=True).fetch(f"{CHROM}:{SIM_S}-{SIM_E}")
        grid[(cell, "experiment")] = hic
        for cond, params in MODELS:
            grid[(cell, cond)] = sim_map(cell, cond, params)
    disp = {k: norm50(crop(m)) for k, m in grid.items() if m is not None}
    allv = np.concatenate([v[np.isfinite(v) & (v > 0)] for v in disp.values()])
    norm = LogNorm(*np.percentile(allv, [2, 99.5]))

    cols = ["experiment"] + [c for c, _ in MODELS]
    nr, nc = len(CELLS), len(cols)
    w = 1.9
    fig = plt.figure(figsize=(w * nc + 1.0, w * nr + 1.6))
    gs = fig.add_gridspec(nr + 1, nc, height_ratios=[1] * nr + [0.28], left=0.06, right=0.91,
                          top=0.9, bottom=0.08, wspace=0.05, hspace=0.07)
    ext = (WS, WE, WE, WS)
    ctr = lambda bp: (bp // BIN) * BIN + BIN / 2
    ticks = np.arange(-(-WS // 500_000) * 500_000, WE + 1, 500_000)
    im = None
    for r, (cell, hic_label, _) in enumerate(CELLS):
        for k, col in enumerate(cols):
            ax = fig.add_subplot(gs[r, k])
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(True); s.set_linewidth(0.4)
            if r == 0:
                ax.set_title("Hi-C" if col == "experiment" else style.MODEL[col]["label"], fontsize=7, pad=4,
                             color=style.TEXT if col == "experiment" else style.MODEL[col]["color"])
            if k == 0:
                ax.set_ylabel(cell, fontsize=9, color=style.CELL[cell], fontweight="bold", labelpad=4)
            m = grid[(cell, col)]
            if m is None:
                ax.set_facecolor("#f3f2ee")
                ax.text(0.5, 0.5, "simulation\npending", ha="center", va="center", fontsize=6.5,
                        color=style.TEXT2, transform=ax.transAxes)
                continue
            im = ax.imshow(disp[(cell, col)], extent=ext, cmap=style.MAP_CMAP, norm=norm,
                           interpolation="none", aspect="equal")
            al.annotate_ep(ax, ctr(al.EBF3_TSS_BP), ctr(al.HS737_BP), "E–P", ext, BIN, oe=ep_oe(m),
                           size_bins=5, fontsize=5.3, text_offset=(-14, -16))
            if col == "experiment":
                ax.text(0.97, 0.03, hic_label, transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=5.5, color=style.TEXT2)

    # gene / enhancer track under every column, same x scale as the maps
    for k in range(nc):
        ax = fig.add_subplot(gs[nr, k])
        ax.set_xlim(WS, WE); ax.set_ylim(-1.2, 1.2)
        ax.set_yticks([])
        for s in ("left", "right", "top"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.4)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t / 1e6:.1f}" for t in ticks], fontsize=5.5)
        ax.tick_params(length=2, pad=1)
        for s, e, strand, name in GENES:
            ebf3 = name == "Ebf3"
            ax.add_patch(plt.Rectangle((max(s, WS), 0.25), min(e, WE) - max(s, WS), 0.45, lw=0,
                                       color=style.TEXT if ebf3 else style.TEXT2))
            if name:
                ax.text((max(s, WS) + min(e, WE)) / 2, 0.95, name, ha="center", va="bottom", fontsize=5.2,
                        color=style.TEXT, fontstyle="italic", fontweight="bold" if ebf3 else "normal")
        ax.plot([al.HS737_BP], [-0.45], marker="^", ms=4, color=style.EP_RING)
        ax.text(al.HS737_BP, -1.15, "hs737", ha="center", va="bottom", fontsize=5.2, color=style.TEXT)
        ax.plot([al.EBF3_TSS_BP], [-0.45], marker="v", ms=4, color=style.EP_RING)
        ax.text(al.EBF3_TSS_BP, -1.15, "TSS", ha="center", va="bottom", fontsize=5.2, color=style.TEXT)
        if k == 0:
            ax.set_xlabel(f"{CHROM} (Mb)", fontsize=6, labelpad=1)

    if im is not None:
        cax = fig.add_axes([0.925, 0.45, 0.012, 0.3])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("contact / mean at 50 kb", fontsize=6.5)
        cb.ax.tick_params(labelsize=5.5)
    fig.suptitle(f"Ebf3 locus, Hi-C vs simulation by cell type — {CHROM}:{WS / 1e6:.2f}–{WE / 1e6:.2f} Mb (10 kb)",
                 y=0.975, fontsize=9.5)
    fig.text(0.06, 0.005,
             "Window: hs737 − 200 kb to Ebf3 TSS + 200 kb. Each map divided by its own mean contact at 50 kb; one shared "
             "log scale. Ring: hs737 × Ebf3-promoter pixel with its O/E (10 kb, diagonal mean over chr7:135.8–137.7 Mb).\n"
             "Models use the row's in vitro CTCF peaks (C, D also its tier-1 enhancer + Ebf3 promoter, 3 kT); "
             "3 replicates × 10,000 conformations pooled. NPC / CN Hi-C is neocortex-sorted (GSE161259 ncx_).",
             fontsize=5.6, color=style.TEXT2, va="bottom")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    style.save(fig, OUT)
    print("wrote", OUT, "and .pdf")


if __name__ == "__main__":
    main()
