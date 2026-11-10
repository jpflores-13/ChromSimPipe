#!/usr/bin/env python
"""
Ebf3 fig1: experimental Hi-C of the simulation window for every available Bonev
2017 state, one shared colour scale, every map annotated with the same function
(scripts/analysis/annotate_loops.py), tracks beneath.

  fig1_experimental_hic.png      full window chr7:135.8-137.7 Mb at 10 kb
  fig1b_experimental_hic_zoom.png  Ebf3 side chr7:137.0-137.7 Mb at 5 kb

States are taken from STATES; in vitro NPC / CN (Stage 1c) are added automatically
when data/mcool/iv_{NPC,CN}_chr7.mcool exist.
Run from the repo root with the ebf3_analysis env.
"""
import gzip, json, os, sys
import numpy as np
import cooler, pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.insert(0, os.getcwd())
from scripts.analysis import style, annotate_loops as al

CHROM = "chr7"
STATES = [("ES", "ES", "ES"), ("ncx_NPC", "NPC", "NPC"), ("ncx_CN", "CN", "CN"),
          ("iv_NPC", "NPC", "NPC"), ("iv_CN", "CN", "CN")]   # (Hi-C, ChIP sample, colour)
STATES = [s for s in STATES if os.path.exists(f"data/mcool/{s[0]}_chr7.mcool")]
TITLE = {"ES": "ES (mESC)", "ncx_NPC": "ncx_NPC (neocortex)", "ncx_CN": "ncx_CN (neocortex)",
         "iv_NPC": "NPC (in vitro)", "iv_CN": "CN (in vitro)"}
EP_RES = 10_000            # pre-registered endpoint resolution (docs/ebf3/01_hic_qc.md)
ATAC_BW = "data/raw/encode_atac/ENCFF097JFT.bigWig"   # ENCODE ENCSR810HQR, E14.5 forebrain: PROXY
style.apply()


def genes(ws, we):
    out = []
    with gzip.open("data/raw/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
        for l in f:
            r = l.split("\t")
            if len(r) > 8 and r[0] == CHROM and r[2] == "gene" and int(r[4]) > ws and int(r[3]) < we \
                    and 'gene_type "protein_coding"' in r[8]:
                out.append((int(r[3]), int(r[4]), r[6], r[8].split('gene_name "')[1].split('"')[0]))
    return out


RNA_STATE = {"ES": "ES", "ncx_NPC": "ncx_NPC", "ncx_CN": "ncx_CN", "iv_NPC": "NPC", "iv_CN": "CN"}
_RNA = {}


def rna_fpkm(state):
    """{gene symbol: mean FPKM of the 2 replicates} from GSE96107 per-GSM tables
    (RNA from the SAME cell source as the Hi-C column)."""
    import glob
    if state not in _RNA:
        reps = []
        for p in sorted(glob.glob(f"data/raw/gse96107_rna/GSM*_{state}_rep*.txt.gz")):
            with gzip.open(p, "rt") as f:
                next(f)
                reps.append({r[2]: float(r[1]) for r in (l.rstrip("\n").split("\t") for l in f) if len(r) > 2})
        _RNA[state] = {g: float(np.mean([d.get(g, np.nan) for d in reps])) for g in reps[0]}
    return _RNA[state]


def ctcf(sample):
    p = f"output_ebf3/loci/ebf3_hs737/ctcf/ctcf_oriented_mm10_{sample}.bed"
    return [((int(r[1]) + int(r[2])) // 2, r[5].strip()) for r in
            (l.split("\t") for l in open(p) if not l.startswith("#"))]


def vista(ws, we):
    import subprocess
    t = subprocess.run(["bigBedToBed", f"-chrom={CHROM}", f"-start={ws}", f"-end={we}",
                        "data/raw/ucsc/vistaEnhancers.bb", "stdout"], capture_output=True, text=True).stdout
    return [(int(r[1]), int(r[2]), r[3]) for r in (l.split("\t") for l in t.splitlines())]


def oe_ep(clr, ws, we, res):
    """E-P O/E with expected = mean of that diagonal inside [ws, we) (same as hic_qc.py)."""
    i, j = al.HS737_BP // res, al.EBF3_TSS_BP // res
    m = clr.matrix(balance=True).fetch(f"{CHROM}:{ws}-{we}")
    i0 = ws // res
    d = j - i
    diag = np.diagonal(m, offset=d)
    return m[i - i0, j - i0] / np.nanmean(diag)


def panel(fig, ax, state, chip, res, ws, we, vmin, vmax, ep_oe, show_cbar=False):
    clr = cooler.Cooler(f"data/mcool/{state}_chr7.mcool::resolutions/{res}")
    m = clr.matrix(balance=True).fetch(f"{CHROM}:{ws}-{we}")
    ext = (ws, we, we, ws)
    im = ax.imshow(m, cmap=style.MAP_CMAP, norm=LogNorm(vmin, vmax), extent=ext,
                   interpolation="none", aspect="equal")
    ax.set_xlim(ws, we); ax.set_ylim(we, ws)
    ax.set_yticks([])
    ax.tick_params(bottom=False, labelbottom=False)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(style.TEXT2); sp.set_linewidth(0.5)
    c = lambda bp: (bp // EP_RES) * EP_RES + EP_RES / 2
    # 5 endpoint bins: at 2-3 bins (50 kb) the ring is invisible on a 1.9 Mb map
    al.annotate_ep(ax, c(al.EBF3_TSS_BP), c(al.HS737_BP), "hs737 × Ebf3", ext, EP_RES, oe=ep_oe,
                   size_bins=5)
    al.mark_ctcf_anchors(ax, ext)
    ax.set_title(TITLE[state], fontsize=8, color=style.TEXT, pad=8)
    div = make_axes_locatable(ax)
    tg = div.append_axes("bottom", size="14%", pad=0.04, sharex=ax)
    te = div.append_axes("bottom", size="12%", pad=0.02, sharex=ax)
    tc = div.append_axes("bottom", size="9%", pad=0.02, sharex=ax)
    tk = div.append_axes("bottom", size="12%", pad=0.02, sharex=ax)
    tx = div.append_axes("bottom", size="12%", pad=0.02, sharex=ax)
    tv = div.append_axes("bottom", size="8%", pad=0.02, sharex=ax)
    ta = div.append_axes("bottom", size="13%", pad=0.02, sharex=ax)
    for a in (tg, te, tc, tk, tx, tv, ta):
        a.set_yticks([]); a.tick_params(labelbottom=False)
        for s in ("left", "right", "top"):
            a.spines[s].set_visible(False)
    span = we - ws
    # fixed rows so neighbouring gene names never collide; 9430038I01Rik unlabelled (caption)
    row = {"Mgmt": 0, "Ebf3": 1, "9430038I01Rik": 0, "Glrx3": 2}
    lab = {"Mgmt", "Ebf3", "Glrx3"}
    for s, e, st, n in genes(ws, we):
        y = 1.0 - 0.45 * row.get(n, 0)
        s0, e0 = max(s, ws), min(e, we)
        tg.add_patch(Rectangle((s0, y - 0.08), e0 - s0, 0.16, color=style.TEXT, lw=0))
        if n in lab:
            left = s0 - span * 0.01 > ws + span * 0.08
            tg.text(s0 - span * 0.01 if left else e0 + span * 0.01, y, n, fontsize=5.5, style="italic",
                    ha="right" if left else "left", va="center", clip_on=True)
    tg.set_ylim(-0.1, 1.2)
    # expression: one bar per protein-coding gene over its body, mean FPKM of 2 replicates
    fp = rna_fpkm(RNA_STATE[state])
    for s, e, st, n in genes(ws, we):
        v = fp.get(n, np.nan)
        if np.isfinite(v):
            s0, e0 = max(s, ws), min(e, we)
            te.add_patch(Rectangle((s0, 0), e0 - s0, np.log10(v + 1), color=style.CELL[chip], lw=0))
            if n == "Ebf3":
                te.text((s0 + e0) / 2, np.log10(v + 1) + 0.15, f"Ebf3 {v:.2g}", fontsize=5.5, ha="center", va="bottom")
    te.set_ylim(0, 3.0)
    for x, st in ctcf(chip):
        dx = span * 0.008 * (1 if st == "+" else -1)
        tc.annotate("", xy=(x + dx, 0.5), xytext=(x - dx * 0.2, 0.5),
                    arrowprops=dict(arrowstyle="-|>,head_width=0.2,head_length=0.35",
                                    color=style.CTCF_STRAND[st], lw=0.8))
    tc.set_ylim(0, 1)
    bw = f"data/raw/gse96107/GSE96107_{chip}_H3K27ac.bw"
    if os.path.exists(bw):
        nb = 400
        v = pyBigWig.open(bw).stats(CHROM, ws, we, nBins=nb, type="mean")
        v = np.array([0 if a is None else a for a in v])
        x = np.linspace(ws, we, nb, endpoint=False) + span / nb / 2
        tk.fill_between(x, 0, v, color=style.CELL[chip], lw=0, step="mid")
        tk.set_ylim(0, max(np.percentile(v, 99.8), 1e-6) * 1.1)
    else:
        tk.text((ws + we) / 2, 0.5, "no H3K27ac in GSE96107", ha="center", va="center",
                fontsize=5.5, color=style.TEXT2)
        tk.set_ylim(0, 1)
    nb = 400
    v = pyBigWig.open(ATAC_BW).stats(CHROM, ws, we, nBins=nb, type="mean")
    v = np.array([0 if a is None else a for a in v])
    xx = np.linspace(ws, we, nb, endpoint=False) + span / nb / 2
    tx.fill_between(xx, 0, v, color=style.TEXT2, lw=0, step="mid")
    tx.set_ylim(0, max(np.percentile(v, 99.8), 1e-6) * 1.1)
    for s, e, n in vista(ws, we):
        hs = n == "hs737"
        tv.add_patch(Rectangle((s, 0.15), max(e - s, span * 0.002), 0.7,
                               color=style.TEXT if hs else "#b9b8b3", lw=0))
    tv.set_ylim(0, 1)
    al.annotate_arc(ta, al.HS737_BP, al.EBF3_TSS_BP, height=0.85)
    ta.set_ylim(0, 1.1)
    ta.set_xlim(ws, we)
    ta.tick_params(labelbottom=True)
    step = 400_000 if span > 1_000_000 else 100_000
    ticks = np.arange(-(-ws // step) * step, we + 1, step)
    ta.set_xticks(ticks); ta.set_xticklabels([f"{t / 1e6:.1f}" for t in ticks], fontsize=6)
    return im, (tg, te, tc, tk, tx, tv, ta)


def figure(res, ws, we, out, title):
    n = len(STATES)
    oe = {}
    mats = []
    for st, _, _ in STATES:
        clr = cooler.Cooler(f"data/mcool/{st}_chr7.mcool::resolutions/{res}")
        mats.append(clr.matrix(balance=True).fetch(f"{CHROM}:{ws}-{we}"))
        ep_clr = cooler.Cooler(f"data/mcool/{st}_chr7.mcool::resolutions/{EP_RES}")
        oe[st] = float(oe_ep(ep_clr, 135_800_000, 137_700_000, EP_RES))   # endpoint O/E: always full window
    allv = np.concatenate([m[np.isfinite(m) & (m > 0)] for m in mats])
    vmin, vmax = np.percentile(allv, 2), np.percentile(allv, 99.5)
    W = 2.55 * n + 0.6
    fig, axes = plt.subplots(1, n, figsize=(W, 6.2))
    fig.subplots_adjust(left=0.08, right=0.9, top=0.86, bottom=0.2, wspace=0.12)
    axes = np.atleast_1d(axes)
    for k, (ax, (st, chip, _)) in enumerate(zip(axes, STATES)):
        im, tr = panel(fig, ax, st, chip, res, ws, we, vmin, vmax, oe[st])
        if k == 0:
            for a, lab in zip(tr, ("genes", "RNA\nlog10 FPKM", "CTCF", "H3K27ac", "ATAC*", "VISTA", "E–P")):
                a.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=6)
    cax = fig.add_axes([0.92, 0.55, 0.012, 0.3])
    cb = fig.colorbar(im, cax=cax); cb.set_label("ICE-balanced contact", fontsize=6); cb.ax.tick_params(labelsize=5.5)
    fig.suptitle(title, y=0.985, fontsize=9)
    fig.legend(handles=al.legend_handles(), loc="upper center", ncol=2, fontsize=6.5, bbox_to_anchor=(0.5, 0.96))
    import textwrap
    cap = (f"Hi-C: GSE161259 (HiCUP reanalysis of Bonev 2017), mm10, {res // 1000} kb, ICE (cis); one shared colour "
             f"scale. Ring: hs737 × Ebf3-promoter pixel at the pre-registered {EP_RES // 1000} kb endpoint resolution, "
             "O/E = pixel / mean of its diagonal within chr7:135.8–137.7 Mb.\n"
             "Tracks: CTCF = GSE96107 IDR peaks oriented by FIMO; H3K27ac = GSE96107 bigWig. ncx_ (neocortex) Hi-C is "
             "paired with IN VITRO NPC / CN ChIP: a cell-source MISMATCH (Ebf3 30–79 FPKM in vitro CN vs 1–4 in ncx_CN). "
             "ES has no H3K27ac in GSE96107. RNA: GSE96107 gene-level FPKM (mean of 2 replicates) from the SAME cell source as each Hi-C column; bar = log10(FPKM+1) over the gene body, Ebf3 value printed.\n*ATAC = PROXY: ENCODE ENCSR810HQR bulk E14.5 forebrain "
             "(fold change, ENCFF097JFT), identical in every column; not from these cells. Triangles: conserved candidate CTCF loop anchors. VISTA: hs737 black, "
             "other elements grey. Unlabelled gene between Ebf3 and Glrx3: 9430038I01Rik. White lines: bins removed by "
             "ICE balancing (low coverage).")
    cap = "\n".join(textwrap.fill(p, width=int(W * 26)) for p in cap.split("\n"))
    fig.text(0.01, 0.005, cap, fontsize=5.8, color=style.TEXT2, va="bottom")
    style.save(fig, out)
    print("wrote", out, {k: round(v, 2) for k, v in oe.items()})


if __name__ == "__main__":
    os.makedirs("figures/ebf3", exist_ok=True)
    figure(10_000, 135_800_000, 137_700_000, "figures/ebf3/fig1_experimental_hic.png",
           "Experimental Hi-C, Ebf3 / hs737 locus, chr7:135.8–137.7 Mb (10 kb)")
    figure(5_000, 137_000_000, 137_700_000, "figures/ebf3/fig1b_experimental_hic_zoom.png",
           "Ebf3 side, chr7:137.0–137.7 Mb (5 kb); the hs737 × Ebf3 pixel lies outside this view")
