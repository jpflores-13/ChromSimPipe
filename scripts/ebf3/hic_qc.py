#!/usr/bin/env python
"""
Ebf3 Stage 1d: Hi-C QC go/no-go gate.

Per cooler: valid pairs (and cis/trans if the pairs file has trans), coverage
across the window, P(s) + local slope, the raw / balanced / O/E value of the
hs737 x Ebf3-promoter pixel at 1, 5, 10, 25 kb, and insulation across
chr7:135.4-138.1 Mb with both candidate window edges marked.

PRE-REGISTERED (written before any count was seen):
  * Endpoint resolution = the finest of 1/5/10/25 kb at which the RAW count in
    the single hs737 x Ebf3-promoter pixel is >= MIN_RAW (30; Poisson CV <= 18%)
    in EVERY state. Fixed once chosen; not changed after simulations are seen.
  * O/E = balanced pixel / mean balanced value of the same diagonal, the mean
    taken INSIDE the simulation window only (the same definition is applied to
    simulated maps, which only exist inside the window).
  * Window-edge check: an edge "sits in a strong boundary" if a cooltools
    insulation boundary with strength above the chr7 75th percentile lies within
    EDGE_TOL (20 kb) of it, at the 10 kb / 200 kb-window insulation track.

Run with the ebf3_analysis env, from the repo root:
  python -I scripts/ebf3/hic_qc.py
"""
import json, os, sys
import numpy as np
import pandas as pd
import cooler, cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())
from scripts.analysis import style

STATES = ["ES", "ncx_NPC", "ncx_CN"]
MCOOL = "data/mcool/{}_chr7.mcool"
CHROM = "chr7"
WIN = (135_800_000, 137_700_000)
WIDE = (135_467_000, 137_930_000)
VIEW = (135_400_000, 138_100_000)
HS737 = 136_083_776          # VISTA hs737 midpoint (element 136,083,204-136,084,349)
TSS = 137_314_445            # Ebf3 TSS (minus strand)
RES_EP = [1000, 5000, 10000, 25000]
MIN_RAW = 30
EDGE_TOL = 20_000
INS_RES, INS_WINDOWS = 10_000, [100_000, 200_000, 500_000]
IGNORE_DIAGS = 2             # cooler balance default; zoomify --balance does not store it in the weights
CAPTION_NOTE = ("ncx_NPC / ncx_CN = GSE161259 reanalysis of neocortex-sorted NPC/CN (GSM2533835-42); "
                "NOT the in vitro NPC/CN that the GSE96107 ChIP and RNA come from.")


def clr(state, res):
    return cooler.Cooler(f"{MCOOL.format(state)}::resolutions/{res}")


def ep_pixel(c, res):
    i, j = HS737 // res, TSS // res
    off = c.offset(CHROM)
    raw = c.matrix(balance=False)[off + i, off + j]
    bal = c.matrix(balance=True)[off + i, off + j]
    lo, hi = WIN[0] // res, WIN[1] // res
    d = j - i
    diag = np.array([c.matrix(balance=True)[off + k, off + k + d] for k in range(lo, hi - d)], float)
    exp_win = np.nanmean(diag)
    return dict(res=res, bin_i=int(i), bin_j=int(j), sep_bins=int(d), raw=int(raw),
                balanced=float(bal), expected_in_window=float(exp_win),
                oe=float(bal / exp_win) if exp_win > 0 else np.nan,
                n_diag_pixels=int(np.isfinite(diag).sum()))


def main():
    style.apply()
    os.makedirs("figures/ebf3", exist_ok=True)
    os.makedirs("results/ebf3", exist_ok=True)
    out = {"pre_registered": dict(min_raw=MIN_RAW, edge_tol_bp=EDGE_TOL, ins_res=INS_RES,
                                  oe_expected="diagonal mean inside chr7:135,800,000-137,700,000"),
           "states": {}}
    ps, cov, ins = {}, {}, {}
    for s in STATES:
        info = {}
        c1 = clr(s, 1000)
        info["resolutions"] = sorted(int(r.split("/")[-1]) for r in cooler.fileops.list_coolers(MCOOL.format(s)))
        info["sum_contacts_1kb_cooler"] = int(c1.info["sum"])
        pairs_stats = f"results/ebf3/pairs_stats_{s}.json"
        if os.path.exists(pairs_stats):
            info.update(json.load(open(pairs_stats)))
        # coverage: raw marginals at 10 kb over the view, and filtered (NaN weight) 1 kb bins in the window
        c10 = clr(s, 10000)
        reg = f"{CHROM}:{VIEW[0]}-{VIEW[1]}"
        m = c10.matrix(balance=False).fetch(reg)
        cov[s] = m.sum(axis=0)
        w1 = c1.bins().fetch(f"{CHROM}:{WIN[0]}-{WIN[1]}")["weight"].values
        info["window_1kb_bins_filtered_frac"] = float(np.isnan(w1).mean())
        info["window_raw_contacts_10kb"] = int(np.triu(c10.matrix(balance=False).fetch(f"{CHROM}:{WIN[0]}-{WIN[1]}")).sum())
        # E-P pixel
        info["ep_pixel"] = [ep_pixel(clr(s, r), r) for r in RES_EP]
        # P(s) on chr7 at 5 kb (balanced), with local slope
        c5 = clr(s, 5000)
        view = pd.DataFrame({"chrom": [CHROM], "start": [0], "end": [c5.chromsizes[CHROM]], "name": [CHROM]})
        cvd = cooltools.expected_cis(c5, view_df=view, smooth=True, aggregate_smoothed=True,
                                     clr_weight_name="weight", ignore_diags=IGNORE_DIAGS, nproc=4)
        cvd["s_bp"] = cvd["dist"] * 5000
        ps[s] = cvd[cvd["dist"] > 0]
        # insulation at 10 kb across chr7
        ins[s] = cooltools.insulation(c10, INS_WINDOWS, view_df=view.assign(end=c10.chromsizes[CHROM]),
                                      clr_weight_name="weight", ignore_diags=IGNORE_DIAGS, nproc=4)
        out["states"][s] = info

    # ---- endpoint resolution (pre-registered rule) ----
    chosen = None
    for k, r in enumerate(RES_EP):
        if all(out["states"][s]["ep_pixel"][k]["raw"] >= MIN_RAW for s in STATES):
            chosen = r
            break
    out["endpoint_resolution_bp"] = chosen

    # ---- window-edge boundary check ----
    edges = {"narrow_left": WIN[0], "narrow_right": WIN[1], "wide_left": WIDE[0], "wide_right": WIDE[1]}
    edge_check = {}
    wkey = "200000"
    for s in STATES:
        t = ins[s]
        bs = t[f"boundary_strength_{wkey}"]
        thr = float(np.nanpercentile(bs, 75))
        strong = t[(t[f"is_boundary_{wkey}"]) & (bs > thr)]
        edge_check[s] = {}
        for name, x in edges.items():
            near = strong[(strong["start"] < x + EDGE_TOL) & (strong["end"] > x - EDGE_TOL)]
            edge_check[s][name] = dict(strong_boundary_within_tol=bool(len(near)),
                                       boundaries=[dict(start=int(a), strength=round(float(b), 3))
                                                   for a, b in zip(near["start"], near[f"boundary_strength_{wkey}"])])
        edge_check[s]["strength_p75"] = round(thr, 3)
    out["edge_check"] = edge_check
    out["widen_window"] = any(edge_check[s][e]["strong_boundary_within_tol"]
                              for s in STATES for e in ("narrow_left", "narrow_right"))

    # ---- table ----
    rows = []
    for s in STATES:
        for p in out["states"][s]["ep_pixel"]:
            rows.append(dict(state=s, **p))
    tab = pd.DataFrame(rows)
    tab.to_csv("results/ebf3/hic_qc_ep_pixel.tsv", sep="\t", index=False)
    with open("results/ebf3/hic_qc.json", "w") as f:
        json.dump(out, f, indent=1, default=float)

    # ---- figure 01: coverage + P(s) + slope + E-P table ----
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.15, 0.8], hspace=0.55, wspace=0.3)
    ax = fig.add_subplot(gs[0, :])
    x = np.arange(VIEW[0], VIEW[1], 10000)[: len(cov[STATES[0]])] / 1e6
    for s in STATES:
        ax.plot(x, cov[s] / np.nanmedian(cov[s]), color=style.CELL[s], lw=1)
        ax.text(x[-1] + 0.02, (cov[s] / np.nanmedian(cov[s]))[-1], s, color=style.TEXT, fontsize=7, va="center")
    for e in WIN:
        ax.axvline(e / 1e6, color=style.TEXT2, lw=0.8)
    for e in WIDE:
        ax.axvline(e / 1e6, color=style.TEXT2, lw=0.8, ls=":")
    for xx, lab in ((HS737, "hs737"), (TSS, "Ebf3 TSS")):
        ax.axvline(xx / 1e6, color=style.EP_RING, lw=0.6, ls="--")
        ax.text(xx / 1e6, ax.get_ylim()[1], lab, fontsize=6.5, ha="center", va="bottom")
    ax.set(xlabel=f"{CHROM} (Mb)", ylabel="raw coverage / median", title="a  Coverage, 10 kb bins (solid: sim window; dotted: widened)")
    ax1, ax2 = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    for s in STATES:
        d = ps[s]
        y = d["balanced.avg.smoothed.agg"] if "balanced.avg.smoothed.agg" in d else d["balanced.avg"]
        y = y / np.interp(1e5, d["s_bp"], y)
        ax1.loglog(d["s_bp"], y, color=style.CELL[s], label=s)
        ls, ly = np.log(d["s_bp"].values), np.log(np.asarray(y, float))
        ax2.semilogx(d["s_bp"].values[1:], np.diff(ly) / np.diff(ls), color=style.CELL[s], label=s)
    for a_ in (ax1, ax2):
        a_.axvline(TSS - HS737, color=style.TEXT2, lw=0.8, ls="--")
        a_.set_xlim(1e4, 1e7)
    ax1.set(xlabel="s (bp)", ylabel="P(s) / P(100 kb)", title="b  P(s), chr7, 5 kb, ICE")
    ax2.set(xlabel="s (bp)", ylabel="d log P / d log s", title="c  Local slope")
    ax1.legend(); ax1.text(TSS - HS737, ax1.get_ylim()[0] * 1.5, " E–P 1.23 Mb", fontsize=6.5, color=style.TEXT2)
    axt = fig.add_subplot(gs[2, :]); axt.axis("off")
    piv = tab.pivot(index="state", columns="res", values="raw").loc[STATES]
    pivoe = tab.pivot(index="state", columns="res", values="oe").loc[STATES]
    cell = [[f"{piv.loc[s, r]:,}  (O/E {pivoe.loc[s, r]:.2f})" for r in RES_EP] for s in STATES]
    tb = axt.table(cellText=cell, rowLabels=STATES, colLabels=[f"{r // 1000} kb" for r in RES_EP],
                   loc="center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(7); tb.scale(1, 1.3)
    axt.set_title(f"d  hs737 × Ebf3-promoter pixel: raw count (O/E, expected = in-window diagonal). "
                  f"Endpoint resolution (rule: raw ≥ {MIN_RAW} in all states): "
                  f"{'%d kb' % (chosen // 1000) if chosen else 'NONE meets the rule'}", fontsize=7.5)
    fig.text(0.01, 0.005, CAPTION_NOTE, fontsize=6, color=style.TEXT2, wrap=True)
    style.save(fig, "figures/ebf3/01_hic_qc.png")

    # ---- figure 02: insulation with both candidate edges ----
    fig, axes = plt.subplots(len(STATES), 1, figsize=(7.2, 5.2), sharex=True)
    for a_, s in zip(axes, STATES):
        t = ins[s]
        t = t[(t["start"] >= VIEW[0]) & (t["end"] <= VIEW[1])]
        xm = (t["start"] + t["end"]) / 2e6
        for w, alpha in zip(INS_WINDOWS, (0.45, 1.0, 0.45)):
            a_.plot(xm, t[f"log2_insulation_score_{w}"], color=style.CELL[s], lw=1 if w == 200000 else 0.7,
                    alpha=alpha)
        b = t[t["is_boundary_200000"] & (t["boundary_strength_200000"] > edge_check[s]["strength_p75"])]
        a_.plot((b["start"] + b["end"]) / 2e6, b["log2_insulation_score_200000"], "v", ms=4,
                color=style.TEXT, label="strong boundary (200 kb, > chr7 p75)")
        for e in WIN:
            a_.axvline(e / 1e6, color=style.TEXT, lw=1)
            a_.axvspan((e - EDGE_TOL) / 1e6, (e + EDGE_TOL) / 1e6, color=style.GRID, zorder=0)
        for e in WIDE:
            a_.axvline(e / 1e6, color=style.TEXT2, lw=1, ls=":")
        a_.set_ylabel(f"{s}\nlog2 insulation")
    axes[0].legend(loc="lower left", fontsize=6.5)
    axes[0].set_title("Insulation (10 kb; 100/200/500 kb windows, 200 kb bold). "
                      "Solid: 135.80 / 137.70 Mb (±20 kb shaded); dotted: widened 135.467 / 137.93 Mb")
    axes[-1].set_xlabel(f"{CHROM} (Mb)")
    fig.text(0.01, -0.01, CAPTION_NOTE, fontsize=6, color=style.TEXT2)
    style.save(fig, "figures/ebf3/02_insulation_window_check.png")
    print(json.dumps({k: out[k] for k in ("endpoint_resolution_bp", "widen_window", "edge_check")}, indent=1))
    print(tab.to_string(index=False))
    print(json.dumps({s: {k: v for k, v in out["states"][s].items() if k != "ep_pixel"} for s in STATES}, indent=1))


if __name__ == "__main__":
    main()
