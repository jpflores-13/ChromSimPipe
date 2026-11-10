#!/usr/bin/env python
"""
Window-edge insulation check for any locus (same rule as the Ebf3 Stage 1d gate):
an edge FAILS if a 200 kb-window insulation boundary (10 kb bins) stronger than
the chromosome's 75th percentile lies within +-EDGE_TOL (20 kb) of it.
Also writes an insulation figure with the edges marked.

  python -I scripts/ebf3/window_edge_check.py --chrom chr2 --start 108800000 --end 110300000 \
      --mcool 'data/mcool/{state}_chr2.mcool' --out results/bdnf/window_edge_check.json \
      --fig figures/bdnf/02_insulation_window_check.png
Exit code 1 if any edge fails in any state.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import cooler, cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd())
from scripts.analysis import style

STATES = ["ES", "ncx_NPC", "ncx_CN"]
EDGE_TOL, RES, W = 20_000, 10_000, 200_000

ap = argparse.ArgumentParser()
ap.add_argument("--chrom", required=True)
ap.add_argument("--start", type=int, required=True)
ap.add_argument("--end", type=int, required=True)
ap.add_argument("--mcool", required=True, help="path pattern with {state}")
ap.add_argument("--out", required=True)
ap.add_argument("--fig", required=True)
a = ap.parse_args()
style.apply()
os.makedirs(os.path.dirname(a.out), exist_ok=True)
os.makedirs(os.path.dirname(a.fig), exist_ok=True)

res, failed = {}, False
fig, axes = plt.subplots(len(STATES), 1, figsize=(7.2, 4.6), sharex=True)
pad = 400_000
for ax, s in zip(axes, STATES):
    c = cooler.Cooler(f"{a.mcool.format(state=s)}::resolutions/{RES}")
    view = pd.DataFrame({"chrom": [a.chrom], "start": [0], "end": [c.chromsizes[a.chrom]], "name": [a.chrom]})
    t = cooltools.insulation(c, [W], view_df=view, clr_weight_name="weight", ignore_diags=2, nproc=4)
    bs = t[f"boundary_strength_{W}"]
    p75 = float(np.nanpercentile(bs, 75))
    strong = t[t[f"is_boundary_{W}"] & (bs > p75)]
    res[s] = {"strength_p75": round(p75, 3)}
    for name, x in (("left", a.start), ("right", a.end)):
        near = strong[(strong["start"] < x + EDGE_TOL) & (strong["end"] > x - EDGE_TOL)]
        nearest = strong.iloc[(strong["start"] + RES / 2 - x).abs().argsort()[:1]]
        res[s][name] = dict(fail=bool(len(near)),
                            nearest_strong_boundary=int(nearest["start"].iloc[0]) if len(nearest) else None,
                            ins_at_edge=round(float(t[t["start"] == x // RES * RES][f"log2_insulation_score_{W}"].iloc[0]), 3))
        failed |= bool(len(near))
    v = t[(t["start"] >= a.start - pad) & (t["end"] <= a.end + pad)]
    ax.plot((v["start"] + RES / 2) / 1e6, v[f"log2_insulation_score_{W}"], color=style.CELL[s], lw=1)
    b = strong[(strong["start"] >= a.start - pad) & (strong["end"] <= a.end + pad)]
    ax.plot((b["start"] + RES / 2) / 1e6, b[f"log2_insulation_score_{W}"], "v", ms=4, color=style.TEXT)
    for x in (a.start, a.end):
        ax.axvline(x / 1e6, color=style.TEXT, lw=1)
        ax.axvspan((x - EDGE_TOL) / 1e6, (x + EDGE_TOL) / 1e6, color=style.GRID, zorder=0)
    ax.set_ylabel(f"{s}\nlog2 insulation")
axes[0].set_title(f"Insulation (10 kb, 200 kb window) around {a.chrom}:{a.start:,}-{a.end:,}; "
                  f"triangles = strong boundaries (> chromosome p75); shaded = ±20 kb of each edge", fontsize=7.5)
axes[-1].set_xlabel(f"{a.chrom} (Mb)")
fig.text(0.01, -0.02, "ncx_ = neocortex Hi-C (GSE161259).", fontsize=6, color=style.TEXT2)
style.save(fig, a.fig)
res["any_edge_fails"] = failed
json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps(res, indent=1))
sys.exit(1 if failed else 0)
