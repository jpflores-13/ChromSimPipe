#!/usr/bin/env python
"""Ebf3 Stage 2 figure: FIMO-oriented CTCF sites per cell type across the
simulation window, as strand-coded arrows; adjacent convergent (+ ... -) pairs
bracketed; hs737 and the Ebf3 TSS marked and linked by an arc.
Run from the repo root: python -I scripts/ebf3/plot_ctcf_orientation.py"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.getcwd())
from scripts.analysis import style, annotate_loops as al

WIN = (135_800_000, 137_700_000)
D = json.load(open("results/ebf3/ctcf_barriers.json"))
style.apply()


def sites(t):
    out = []
    for l in open(f"output_ebf3/loci/ebf3_hs737/ctcf/ctcf_oriented_mm10_{t}.bed"):
        if l.startswith("#"):
            continue
        r = l.split("\t")
        out.append(((int(r[1]) + int(r[2])) // 2, r[5].strip()))
    return sorted(out)


fig, axes = plt.subplots(4, 1, figsize=(7.2, 3.9), sharex=True,
                         gridspec_kw=dict(height_ratios=[1, 1, 1, 0.9], hspace=0.15))
L = 9_000          # arrow length (bp)
for ax, t in zip(axes, ("ES", "NPC", "CN")):
    s = sites(t)
    for x, st in s:
        dx = L if st == "+" else -L
        ax.annotate("", xy=(x + dx, 0), xytext=(x - dx * 0.2, 0),
                    arrowprops=dict(arrowstyle="-|>,head_width=0.25,head_length=0.5",
                                    color=style.CTCF_STRAND[st], lw=1.2))
    for (a, sa), (b, sb) in zip(s, s[1:]):
        if sa == "+" and sb == "-":
            ax.plot([a, a, b, b], [0.35, 0.6, 0.6, 0.35], color=style.TEXT2, lw=0.8)
    d = D[t]
    ax.text(1.005, 0.5, f"{t}\n{d['n']} sites ({d['plus']}+ / {d['minus']}−)\n"
                        f"{d['convergent_adjacent']} adjacent conv.", transform=ax.transAxes,
            fontsize=6.5, va="center", color=style.TEXT)
    ax.set_ylim(-0.6, 0.9); ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    for xx in (al.HS737_BP, al.EBF3_TSS_BP):
        ax.axvline(xx, color=style.EP_RING, lw=0.6, ls="--", zorder=0)
tr = axes[3]
al.annotate_arc(tr, al.HS737_BP, al.EBF3_TSS_BP, height=0.8, label="hs737 ↔ Ebf3 TSS, 1.23 Mb")
for xx, lab in ((al.HS737_BP, "hs737"), (al.EBF3_TSS_BP, "Ebf3 TSS")):
    tr.plot([xx], [0], marker="o", ms=4, color=style.EP_RING)
    tr.text(xx, -0.3, lab, ha="center", va="top", fontsize=6.5)
tr.plot([a * 1000 for a in al.CTCF_ANCHORS_KB], [0.0] * len(al.CTCF_ANCHORS_KB), "v", ms=4,
        color=style.CTCF_ANCHOR)
tr.set_ylim(-0.9, 1.25); tr.set_yticks([]); tr.spines["left"].set_visible(False)
tr.set_xlim(*WIN)
tr.set_xticks(np.arange(135.8e6, 137.71e6, 0.2e6))
tr.set_xticklabels([f"{v / 1e6:.1f}" for v in tr.get_xticks()])
tr.set_xlabel("chr7 (Mb, mm10)")
axes[0].legend(handles=[Line2D([], [], color=style.CTCF_STRAND["+"], marker=">", lw=1, label="+ motif (stops leftward legs)"),
                        Line2D([], [], color=style.CTCF_STRAND["-"], marker="<", lw=1, label="− motif (stops rightward legs)"),
                        Line2D([], [], color=style.TEXT2, lw=0.8, label="adjacent convergent pair"),
                        al.legend_handles()[1]],
               ncol=4, loc="lower left", bbox_to_anchor=(0, 1.02), fontsize=6.5)
fig.text(0.01, -0.06,
         "CTCF: GSE96107 IDR-0.05 peaks (in vitro ES / NPC / CN), oriented by FIMO MA0139.1 (p < 1e-4, ±250 bp); "
         "peaks without a hit dropped (ES 15, NPC 4, CN 1).\nOnly ADJACENT convergent pairs are bracketed "
         "(all + … − pairs: ES 610, NPC 43, CN 16). CN calls 4.3× fewer CTCF peaks genome-wide than ES "
         "(peak-call sensitivity, not shown to be CTCF loss).", fontsize=6, color=style.TEXT2)
os.makedirs("figures/ebf3", exist_ok=True)
style.save(fig, "figures/ebf3/03_ctcf_orientation.png")
print("wrote figures/ebf3/03_ctcf_orientation.png")
