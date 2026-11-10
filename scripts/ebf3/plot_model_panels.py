#!/usr/bin/env python
"""
Ebf3 fig2 / fig3 and the primary endpoint.

Reads each condition's per-replicate simulated map
  output_ebf3/loci/ebf3_hs737/sims/merged_<params_name>_ctcf-CN_rep<N>/analysis/sim_contact_map.npy
(1 kb, averaged over the 28 tiled copies; written by analyzeCondition).

Primary endpoint (PRE-REGISTERED in scripts/ebf3/hic_qc.py before any simulation
existed): E-P O/E = the single 10 kb hs737 x Ebf3-promoter pixel divided by the
mean of its diagonal inside chr7:135,800,000-137,700,000. Applied identically
to simulated (1 kb maps summed into 10 kb bins) and experimental maps.
CI: bootstrap over replicates (resample replicates with replacement; with 3
replicates the CI is coarse; per-replicate values are also reported).

  figures/ebf3/fig2_model_panels.png   ncx_CN Hi-C + A, A0, B, C, D
  figures/ebf3/fig3_ep_endpoint.png    E-P O/E per model, experimental references
  results/ebf3/ep_endpoint.tsv         per replicate
  results/ebf3/model_comparison.csv    rows metric=ep_oe_10kb (appended / replaced)
Run from the repo root with the ebf3_analysis env.
"""
import glob, json, os, sys
import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
sys.path.insert(0, os.getcwd())
from scripts.analysis import style, annotate_loops as al

WS, WE, RES, EP_RES = 135_800_000, 137_700_000, 1000, 10_000
F = EP_RES // RES
SIMS = "output_ebf3/loci/ebf3_hs737/sims"
MODELS = [("A_polymer", "polymer"), ("A0_polymer_attr", "polymer_attr"), ("B_extrusion", "extrusion"),
          ("C_sticky_t1", "sticky_t1"), ("D_combined_t1", "combined_t1")]
EXP = "ncx_CN"
NBOOT = 10_000
style.apply()


def coarsen(m, f=F):
    n = m.shape[0] // f
    m = np.nan_to_num(m[: n * f, : n * f])
    return m.reshape(n, f, n, f).sum(axis=(1, 3))


def oe_map(m):
    out = np.full_like(m, np.nan, dtype=float)
    for d in range(m.shape[0]):
        diag = np.diagonal(m, d)
        mu = np.nanmean(diag)
        if mu > 0:
            idx = np.arange(m.shape[0] - d)
            out[idx, idx + d] = diag / mu
            out[idx + d, idx] = diag / mu
    return out


def ep_oe(m10):
    i, j = al.HS737_BP // EP_RES - WS // EP_RES, al.EBF3_TSS_BP // EP_RES - WS // EP_RES
    return float(m10[i, j] / np.nanmean(np.diagonal(m10, j - i)))


def norm50(m10):
    return m10 / np.nanmean(np.diagonal(m10, 50_000 // EP_RES))


def load_models():
    out = {}
    for name, pn in MODELS:
        reps = sorted(glob.glob(f"{SIMS}/merged_{pn}_ctcf-CN_rep*/analysis/sim_contact_map.npy"))
        if reps:
            out[name] = [coarsen(np.load(p)) for p in reps]
    return out


def boot_ci(v, rng):
    v = np.asarray(v, float)
    if len(v) < 2:
        return (np.nan, np.nan)
    b = rng.choice(v, size=(NBOOT, len(v)), replace=True).mean(axis=1)
    return tuple(np.percentile(b, [2.5, 97.5]))


def main():
    rng = np.random.default_rng(20261110)
    models = load_models()
    if not models:
        sys.exit(f"no sim_contact_map.npy under {SIMS}/merged_*/analysis yet")
    exp = {}
    for s in ("ES", "ncx_NPC", "ncx_CN"):
        c = cooler.Cooler(f"data/mcool/{s}_chr7.mcool::resolutions/{EP_RES}")
        exp[s] = c.matrix(balance=True).fetch(f"chr7:{WS}-{WE}")
    halves = json.load(open("results/ebf3/ep_oe_experimental_halves.json"))

    # ---- endpoint table ----
    rows = []
    for name, reps in models.items():
        for k, m in enumerate(reps):
            rows.append(dict(model=name, replicate=k, ep_oe_10kb=ep_oe(m)))
    ep = pd.DataFrame(rows)
    ep.to_csv("results/ebf3/ep_endpoint.tsv", sep="\t", index=False)
    comp = []
    for name, g in ep.groupby("model", sort=False):
        lo, hi = boot_ci(g["ep_oe_10kb"], rng)
        comp.append(dict(model=name, metric="ep_oe_10kb", state="sim (CN CTCF)", value=g["ep_oe_10kb"].mean(),
                         ci_low=lo, ci_high=hi, n_replicates=len(g), kind="prediction",
                         note="bootstrap over replicates"))
    for s, h in halves.items():
        comp.append(dict(model="experiment", metric="ep_oe_10kb", state=s, value=h["full"],
                         ci_low=min(h["half_A"], h["half_B"]), ci_high=max(h["half_A"], h["half_B"]),
                         n_replicates=2, kind="observation", note="range = split-half values (replicate ceiling)"))
    comp = pd.DataFrame(comp)
    path = "results/ebf3/model_comparison.csv"
    if os.path.exists(path):
        old = pd.read_csv(path)
        comp = pd.concat([old[old["metric"] != "ep_oe_10kb"], comp], ignore_index=True)
    comp.to_csv(path, index=False)
    print(comp[comp.metric == "ep_oe_10kb"].round(3).to_string(index=False))

    # ---- fig2 ----
    cols = [("experiment", [exp[EXP]])] + [(n, models[n]) for n, _ in MODELS if n in models]
    disp = {n: norm50(np.mean(r, axis=0)) for n, r in cols}
    oes = {n: oe_map(np.mean(r, axis=0)) for n, r in cols}
    allv = np.concatenate([v[np.isfinite(v) & (v > 0)] for v in disp.values()])
    vmin, vmax = np.percentile(allv, 2), np.percentile(allv, 99.5)
    n = len(cols)
    fig, axes = plt.subplots(2, n, figsize=(1.75 * n + 0.8, 4.6))
    fig.subplots_adjust(left=0.03, right=0.9, top=0.8, bottom=0.18, wspace=0.06, hspace=0.12)
    ext = (WS, WE, WE, WS)
    c = lambda bp: (bp // EP_RES) * EP_RES + EP_RES / 2
    for k, (name, reps) in enumerate(cols):
        val = ep_oe(np.mean(reps, axis=0))
        title = f"Hi-C {EXP}" if name == "experiment" else style.MODEL[name]["label"]
        for row, (data, kw) in enumerate(((disp[name], dict(cmap=style.MAP_CMAP, norm=LogNorm(vmin, vmax))),
                                           (np.log2(oes[name]), dict(cmap=style.OE_CMAP,
                                                                     norm=TwoSlopeNorm(0, -1.5, 1.5))))):
            ax = axes[row, k]
            im = ax.imshow(data, extent=ext, interpolation="none", aspect="equal", **kw)
            ax.set_xticks([]); ax.set_yticks([])
            al.annotate_ep(ax, c(al.EBF3_TSS_BP), c(al.HS737_BP), "E–P", ext, EP_RES, oe=val,
                           size_bins=5, fontsize=5.5, text_offset=(-14, -16))
            al.mark_ctcf_anchors(ax, ext, ms=2.5)
            if row == 0:
                ax.set_title(title, fontsize=6.8, pad=6,
                             color=style.TEXT if name == "experiment" else style.MODEL[name]["color"])
            if k == n - 1:
                cax = fig.add_axes([0.91, 0.52 if row == 0 else 0.2, 0.01, 0.25])
                cb = fig.colorbar(im, cax=cax)
                cb.set_label("contact / mean at 50 kb" if row == 0 else "log2 O/E", fontsize=6)
                cb.ax.tick_params(labelsize=5.5)
    fig.suptitle("Simulated vs experimental contacts, chr7:135.8–137.7 Mb (10 kb)", y=0.985, fontsize=9)
    fig.legend(handles=al.legend_handles(), loc="upper center", ncol=2, fontsize=6.3, bbox_to_anchor=(0.47, 0.94))
    fig.text(0.01, 0.005,
             "Top: each map divided by its own mean contact at 50 kb (sim and Hi-C units differ), one shared log scale. "
             "Bottom: observed/expected (each diagonal / its mean), one shared scale. Models: CN CTCF (in vitro), 3 replicates "
             "pooled, 28 tiled copies each.\nHi-C: GSE161259 ncx_CN (NEOCORTEX), mismatched to the in vitro CN CTCF used by "
             "the models. Compare B with A; C and D with A0 (A0/C/D use polysim's heteropolymer force).",
             fontsize=5.6, color=style.TEXT2, va="bottom")
    style.save(fig, "figures/ebf3/fig2_model_panels.png")

    # ---- fig3 ----
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    xs = {name: i for i, (name, _) in enumerate(MODELS)}
    for s, h in halves.items():
        ax.axhspan(min(h["half_A"], h["half_B"]), max(h["half_A"], h["half_B"]), color=style.CELL[s], alpha=0.18, lw=0)
        ax.axhline(h["full"], color=style.CELL[s], lw=1.2)
        ax.text(len(MODELS) - 0.45, h["full"], f"{s} Hi-C {h['full']:.2f}", fontsize=6.3, va="center", color=style.TEXT)
    for _, r in comp[(comp.metric == "ep_oe_10kb") & (comp.model != "experiment")].iterrows():
        x = xs[r["model"]]
        v = ep[ep.model == r["model"]]["ep_oe_10kb"].values
        ax.plot(x + np.linspace(-0.08, 0.08, len(v)), v, "o", ms=3, mfc="none", mec=style.TEXT2)
        ax.errorbar(x, r["value"], yerr=[[r["value"] - r["ci_low"]], [r["ci_high"] - r["value"]]],
                    fmt="o", ms=5, color=style.TEXT, capsize=2, lw=1)
    ax.axhline(1, color=style.TEXT2, lw=0.6, ls=":")
    ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels([style.MODEL[m]["label"].replace("  ", "\n", 1) for m, _ in MODELS], fontsize=6.3)
    ax.set_xlim(-0.5, len(MODELS) + 0.6)
    ax.set_ylabel("E–P observed / expected (10 kb)")
    ax.set_title("hs737 × Ebf3 promoter, 1.23 Mb", fontsize=8.5)
    fig.text(0.01, -0.12,
             "Black: model mean ± bootstrap 95% CI over replicates (open circles: replicates); predictions with untuned "
             "cohesin defaults. Lines: experimental O/E; shaded: range of the two random halves of each dataset "
             "(replicate ceiling). ncx_ Hi-C is neocortex, models use in vitro CN CTCF.",
             fontsize=5.6, color=style.TEXT2, wrap=True)
    style.save(fig, "figures/ebf3/fig3_ep_endpoint.png")
    print("wrote fig2, fig3")


if __name__ == "__main__":
    main()
