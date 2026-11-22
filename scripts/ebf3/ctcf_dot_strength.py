#!/usr/bin/env python
"""
Per-CTCF-site Hi-C dot strength and manual stall-probability tiers (Michele: 0.35 for
weaker sites, 0.5 for sites with a clear Hi-C dot, 0.6 for exceptionally strong ones).

  python scripts/ebf3/ctcf_dot_strength.py --locus ebf3 --chrom chr7 --start 134830000 --end 138590000 \
      --beds output_ebf3_wide/loci/ebf3_wide/ctcf/ctcf_oriented_mm10_{cell}.bed \
      --hic data/mcool/{hic}_chr7.mcool
  (defaults for --locus ebf3 / bdnf are built in; rerun Bdnf with --start/--end/--beds for a wider region)

Method (10 kb; 5 kb has 1-4 raw counts per pixel at 1 Mb in these data, too noisy):
  * O/E map: balanced (ICE `weight`) counts / cooltools expected_cis over the whole chromosome.
  * Candidate dots: every convergent pair, forward (+) site upstream of a reverse (-) site,
    50 kb - 1.5 Mb apart (bin distance).
  * Enrichment of a pair: mean O/E in the 3x3 peak around the pixel (+-10 kb slop) over the
    mean O/E in a HiCCUPS-style donut (11x11 window minus the 3x3 peak and minus the row and
    column through the pixel). Distance decay is removed by O/E, so this is a local contrast.
    `oe` = mean O/E of the 3x3 peak alone (no local background).
  * Null: each pair moved as a whole by a random +-50-200 kb (same distance, anchors off CTCF),
    100 draws per pair. A site's null is the best enrichment over the same set of partner
    distances in one draw (the max over k partners is biased upward, so the site is compared with
    the max of k null pairs); its clear threshold is the 95th percentile of that.
  * Per site: best-enrichment convergent partner, its distance, and the number of partners
    with enrichment > --partner-min (1.5).
Tiers (per cell type): tier_s = 0.5 if best_enrichment >= the site's own null threshold
(or --clear as one number for every site); 0.6 if also in the top --exceptional-quantile
(default 0.9, i.e. top 10%) of the clear sites' best_enrichment; else 0.35. Sites that share
a bin with a stronger site inherit the bin's statistics (dots are resolved per 10 kb bin).

Writes results/<locus>/ctcf_site_tiers_<cell>.tsv (one row per oriented-BED row) and
figures/<locus>/ctcf_dot_strength.{png,pdf}. Run from the repo root with the ebf3_analysis env.
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import bioframe
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd())
from scripts.analysis import style

CELLS = ("ES", "NPC", "CN")
HIC_NAME = {"ES": "ES", "NPC": "ncx_NPC", "CN": "ncx_CN"}
DEFAULTS = {
    "ebf3": dict(chrom="chr7", start=134_830_000, end=138_590_000,
                 beds="output_ebf3_wide/loci/ebf3_wide/ctcf/ctcf_oriented_mm10_{cell}.bed",
                 hic="data/mcool/{hic}_chr7.mcool"),
    "bdnf": dict(chrom="chr2", start=108_800_000, end=110_300_000,
                 beds="output_bdnf/loci/bdnf/ctcf/ctcf_oriented_mm10_{cell}.bed",
                 hic="data/mcool/{hic}_chr2.mcool"),
}
TIER_COLOR = {0.35: "#b8b6af", 0.5: "#2a78d6", 0.6: "#c0392b"}


def args_():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--locus", required=True)
    ap.add_argument("--chrom"); ap.add_argument("--start", type=int); ap.add_argument("--end", type=int)
    ap.add_argument("--beds", help="oriented CTCF BED pattern with {cell}")
    ap.add_argument("--hic", help="mcool pattern with {hic} (ES, ncx_NPC, ncx_CN)")
    ap.add_argument("--res", type=int, default=10_000)
    ap.add_argument("--min-dist", type=int, default=50_000)
    ap.add_argument("--max-dist", type=int, default=1_500_000)
    ap.add_argument("--peak", type=int, default=1, help="peak half-width in bins (3x3)")
    ap.add_argument("--window", type=int, default=5, help="donut half-width in bins (11x11)")
    ap.add_argument("--partner-min", type=float, default=1.5)
    ap.add_argument("--clear", default="site",
                    help="clear-dot threshold: 'site' = best_enrichment above the site's own null "
                         "(--null-pct of the best-of-same-partner-distances at random positions), or a number")
    ap.add_argument("--null-pct", type=float, default=95.0)
    ap.add_argument("--exceptional-quantile", type=float, default=0.9)
    ap.add_argument("--null-draws", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20261112)
    a = ap.parse_args()
    for k, v in DEFAULTS.get(a.locus, {}).items():
        if getattr(a, k) is None:
            setattr(a, k, v)
    missing = [k for k in ("chrom", "start", "end", "beds", "hic") if getattr(a, k) is None]
    if missing:
        ap.error(f"--locus {a.locus} has no defaults; give {missing}")
    return a


def read_sites(path, chrom, start, end):
    df = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=range(6),
                     names=["chrom", "start", "end", "name", "score", "strand"])
    df = df[(df.chrom == chrom) & (df.end > start) & (df.start < end)].reset_index(drop=True)
    return df


def oe_region(clr, chrom, r0, r1, exp):
    m = clr.matrix(balance=True).fetch(f"{chrom}:{r0}-{r1}")
    n = m.shape[0]
    d = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    e = exp.reindex(range(n))["balanced.avg"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        oe = m / e[d]
    oe[~np.isfinite(oe)] = np.nan
    return oe


def enrichment(oe, i, j, p, w):
    """(local enrichment, peak O/E) for pixel (i, j), i < j, in the O/E map."""
    n = oe.shape[0]
    if i - w < 0 or j + w >= n:
        return np.nan, np.nan
    win = oe[i - w:i + w + 1, j - w:j + w + 1].copy()
    peak = win[w - p:w + p + 1, w - p:w + p + 1]
    pk = np.nanmean(peak) if np.isfinite(peak).sum() >= (2 * p + 1) ** 2 // 2 + 1 else np.nan
    donut = win.copy()
    donut[w - p:w + p + 1, w - p:w + p + 1] = np.nan
    donut[w, :] = np.nan
    donut[:, w] = np.nan
    # pixels too close to the diagonal (near-diagonal O/E is unreliable)
    ii, jj = np.meshgrid(np.arange(i - w, i + w + 1), np.arange(j - w, j + w + 1), indexing="ij")
    donut[(jj - ii) < 3] = np.nan
    if np.isfinite(donut).sum() < 10:
        return np.nan, pk
    bg = np.nanmean(donut)
    return (pk / bg if bg > 0 else np.nan), pk


def analyse_cell(a, cell, rng):
    sites = read_sites(a.beds.format(cell=cell), a.chrom, a.start, a.end)
    clr = cooler.Cooler(f"{a.hic.format(hic=HIC_NAME[cell])}::resolutions/{a.res}")
    view = pd.DataFrame(dict(chrom=[a.chrom], start=[0], end=[clr.chromsizes[a.chrom]], name=[a.chrom]))
    exp = cooltools.expected_cis(clr, view_df=view, smooth=False, ignore_diags=2,
                                 clr_weight_name="weight", nproc=1)
    exp = exp[exp.region1 == a.chrom].set_index("dist")
    pad = (a.window + 25) * a.res          # donut + null shifts beyond the outermost sites
    r0 = max(0, (a.start - pad) // a.res * a.res)
    r1 = min(clr.chromsizes[a.chrom], -(-(a.end + pad) // a.res) * a.res)
    oe = oe_region(clr, a.chrom, r0, r1, exp)
    sites["bin"] = ((sites.start + sites.end) // 2 - r0) // a.res
    fwd, rev = sites[sites.strand == "+"], sites[sites.strand == "-"]
    dmin, dmax = a.min_dist // a.res, a.max_dist // a.res

    pairs = []
    for fi, f in fwd.iterrows():
        for ri, r in rev.iterrows():
            d = r.bin - f.bin
            if dmin <= d <= dmax:
                en, pk = enrichment(oe, int(f.bin), int(r.bin), a.peak, a.window)
                pairs.append(dict(fwd=fi, rev=ri, dist_bp=int(d * a.res), enrichment=en, oe=pk))
    pairs = pd.DataFrame(pairs, columns=["fwd", "rev", "dist_bp", "enrichment", "oe"])

    # null: each pair moved as a whole by +-5..20 bins (same distance, anchors off CTCF);
    # nullmat[pair, draw]
    nullmat = np.full((len(pairs), a.null_draws), np.nan)
    for k, pr in enumerate(pairs.itertuples()):
        i = int(sites.loc[int(pr.fwd), "bin"]); d = int(pr.dist_bp) // a.res
        for t in range(a.null_draws):
            s = int(rng.integers(5, 21) * rng.choice([-1, 1]))
            nullmat[k, t] = enrichment(oe, i + s, i + s + d, a.peak, a.window)[0]
    null = nullmat[np.isfinite(nullmat)]

    rows = []
    for si, s in sites.iterrows():
        col, other = ("fwd", "rev") if s.strand == "+" else ("rev", "fwd")
        sub = pairs[(pairs[col] == si) & np.isfinite(pairs.enrichment)]
        # per-site null: the best of the same partner distances at random positions, per draw
        nmax = np.nanmax(nullmat[sub.index.to_numpy()], axis=0) if len(sub) else np.array([np.nan])
        site_null = float(np.nanpercentile(nmax, a.null_pct)) if np.isfinite(nmax).any() else np.nan
        if len(sub):
            b = sub.loc[sub.enrichment.idxmax()]
            rows.append(dict(best_enrichment=float(b.enrichment), best_oe=float(b.oe),
                             best_partner=sites.loc[int(b[other]), "name"], partner_dist_bp=int(b.dist_bp),
                             n_partners=int((sub.enrichment > a.partner_min).sum()), site_null=site_null))
        else:
            rows.append(dict(best_enrichment=np.nan, best_oe=np.nan, best_partner="NA",
                             partner_dist_bp=np.nan, n_partners=0, site_null=np.nan))
    out = pd.concat([sites[["name", "chrom", "start", "end", "strand"]], pd.DataFrame(rows)], axis=1)

    if str(a.clear) == "site":           # per-site null (default)
        is_clear = out.best_enrichment >= out.site_null
        clear = float(np.nanmedian(out.site_null))          # summary only
    else:
        clear = float(a.clear)
        is_clear = out.best_enrichment >= clear
    exc = (float(np.quantile(out.best_enrichment[is_clear], a.exceptional_quantile))
           if is_clear.sum() >= 2 else np.inf)
    tier = np.where(is_clear, 0.5, 0.35)
    tier = np.where(is_clear & (out.best_enrichment >= exc) & (is_clear.sum() >= 2), 0.6, tier)
    out["tier_s"] = tier
    return out, pairs, null, dict(clear=clear, exceptional=exc, n_pairs=len(pairs), n_null=len(null))


def main():
    a = args_()
    rng = np.random.default_rng(a.seed)
    style.apply()
    os.makedirs(f"results/{a.locus}", exist_ok=True)
    os.makedirs(f"figures/{a.locus}", exist_ok=True)
    res, thresholds = {}, {}
    for cell in CELLS:
        out, pairs, null, th = analyse_cell(a, cell, rng)
        cols = ["name", "chrom", "start", "end", "strand", "best_enrichment", "best_partner",
                "partner_dist_bp", "n_partners", "tier_s", "best_oe", "site_null"]
        out[cols].to_csv(f"results/{a.locus}/ctcf_site_tiers_{cell}.tsv", sep="\t", index=False,
                         float_format="%.4g")
        res[cell] = (out, pairs, null)
        thresholds[cell] = th
        print(f"{a.locus} {cell}: {len(out)} sites, {th['n_pairs']} convergent pairs; "
              f"clear >= {th['clear']:.3f}, exceptional >= {th['exceptional']:.3f}; tiers "
              + ", ".join(f"{t}: {int((out.tier_s == t).sum())}" for t in (0.35, 0.5, 0.6)))
    json.dump(dict(params={k: v for k, v in vars(a).items()}, thresholds=thresholds),
              open(f"results/{a.locus}/ctcf_site_tiers_params.json", "w"), indent=1, default=str)

    # ---- QC figure: enrichment distributions (top), sites along the locus (bottom) ----
    fig, axes = plt.subplots(2, 3, figsize=(11, 5.6), gridspec_kw=dict(height_ratios=[1, 1.1]))
    for k, cell in enumerate(CELLS):
        out, pairs, null = res[cell]; th = thresholds[cell]
        ax = axes[0, k]
        bins = np.linspace(0, max(4, np.nanmax(np.r_[pairs.enrichment, null, 1])), 50)
        ax.hist(null, bins=bins, density=True, color="#d6d5d0", label="null (pair moved off CTCF)")
        ax.hist(pairs.enrichment.dropna(), bins=bins, density=True, histtype="step", color=style.TEXT,
                lw=1, label="convergent CTCF pairs")
        ax.hist(out.best_enrichment.dropna(), bins=bins, density=True, histtype="step",
                color=style.CELL[cell], lw=1.6, label="per-site best")
        ax.axvline(th["clear"], color=TIER_COLOR[0.5], ls="--", lw=1)   # median per-site threshold
        if np.isfinite(th["exceptional"]):
            ax.axvline(th["exceptional"], color=TIER_COLOR[0.6], ls="--", lw=1)
        ax.set_title(f"{cell}  (Hi-C {HIC_NAME[cell]})", color=style.CELL[cell], fontsize=9)
        ax.set_xlabel("local enrichment (peak / donut, O/E)")
        if k == 0:
            ax.set_ylabel("density"); ax.legend(fontsize=6)
        ax = axes[1, k]
        for t, g in out.groupby("tier_s"):
            for strand, mk in (("+", ">"), ("-", "<")):
                h = g[g.strand == strand]
                ax.scatter((h.start + h.end) / 2e6, h.best_enrichment.fillna(0), marker=mk, s=22,
                           color=TIER_COLOR[t], label=f"s = {t}" if strand == "+" else None, lw=0)
        ax.scatter((out.start + out.end) / 2e6, out.site_null, marker="_", s=30, color=TIER_COLOR[0.5],
                   lw=0.8, label="site null (95%)" if k == 0 else None)
        ax.set_xlim(a.start / 1e6, a.end / 1e6)
        ax.set_xlabel(f"{a.chrom} (Mb)")
        if k == 0:
            ax.set_ylabel("best convergent-partner enrichment"); ax.legend(fontsize=6)
    fig.suptitle(f"{a.locus}: CTCF-site Hi-C dot strength ({a.res // 1000} kb) and stall tiers", fontsize=10)
    fig.tight_layout()
    style.save(fig, f"figures/{a.locus}/ctcf_dot_strength.png")


if __name__ == "__main__":
    main()
