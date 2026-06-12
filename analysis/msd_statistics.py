#!/usr/bin/env python
"""
Between-condition statistics for two-point MSD fits.

What this module does
---------------------
Each replicate's MSD analysis produces one scalar for each fitted
parameter per probe pair: the anomalous exponent ``alpha`` and the
generalised diffusion coefficient ``K_alpha`` (so MSD(tau) = K_alpha *
tau^alpha).  With ``n_rep`` replicates per condition, we obtain an
``n_rep``-long vector of alphas and K_alphas per (condition, pair).

This module compares two conditions on those vectors and reports a
panel of diagnostics:

  - mean and SEM of each condition
  - the difference of means
  - 95 percent bootstrap confidence interval for the difference
  - Welch's two-sided t-test p-value and degrees of freedom
  - Mann-Whitney U two-sided p-value (rank-based, distribution-free)
  - a label-permutation two-sided p-value (exact when feasible)
  - Cohen's d effect size using the pooled within-condition std

Why all four tests instead of just one
--------------------------------------
With ``n_rep`` = 3 per condition, no single frequentist test carries
much weight: Welch's t assumes approximate normality, Mann-Whitney is
valid but extremely underpowered at this sample size, and the
permutation test becomes exact but has only C(6,3) = 20 possible
label arrangements so the smallest achievable two-sided p is ~0.1.

Report them together, because their agreement (or disagreement) is
more informative than any single p-value.  When the three disagree in
a small-sample setting, the effect size and the bootstrap CI are the
numbers that should guide interpretation: they quantify *how much* the
conditions differ and *how precisely* you have measured that
difference.  A p-value alone cannot substitute for either.

What NOT to do
--------------
Do not treat per-tile fits from the same replicate as independent
replicates for a between-condition test.  Tiles share one polymer
chain and one LEF pool; they are pseudo-replicates, not independent
samples.  They are useful for within-replicate uncertainty bands, and
the standalone per-tile fit is stored alongside the pooled fit for
exactly that purpose, but the comparison functions in this module
deliberately operate on per-replicate estimates only.

Public API
----------
compare_alpha          : one-shot pairwise comparison on alpha vectors.
compare_K_alpha        : same, on K_alpha (log-transformed internally).
collect_per_replicate  : aggregate per-replicate JSONs for one pair.
run_msd_statistics     : orchestrator that walks a common_dir of
                         pooled outputs, finds every (pair, condition)
                         combination, and emits a summary JSON and a
                         forest plot.
plot_forest            : forest plot (mean ± SE or CI) of one parameter.
"""

from __future__ import annotations

import glob
import itertools
import json
import logging
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# CORE PAIRWISE STATISTICS
# =============================================================================

def _bootstrap_diff_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int = 10_000,
    q_lo: float = 2.5,
    q_hi: float = 97.5,
    seed: Optional[int] = 0,
) -> Tuple[float, float, float]:
    """
    Percentile bootstrap CI for the difference of means ``mean(x) - mean(y)``.
    Resamples with replacement within each group and recomputes the
    difference ``n_boot`` times.  Returns (diff, ci_lo, ci_hi).

    Works even for ``n=3`` per group but its variance is high; interpret
    the width of the CI as a rough scale of uncertainty rather than a
    strictly calibrated coverage statement.
    """
    rng = np.random.default_rng(seed)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return float("nan"), float("nan"), float("nan")
    ix = rng.integers(0, nx, size=(n_boot, nx))
    iy = rng.integers(0, ny, size=(n_boot, ny))
    diffs = x[ix].mean(axis=1) - y[iy].mean(axis=1)
    return (float(np.mean(x) - np.mean(y)),
            float(np.percentile(diffs, q_lo)),
            float(np.percentile(diffs, q_hi)))


def _welch_t(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """
    Welch's two-sample t-test (two-sided).  Uses scipy if available,
    otherwise computes the statistic and a normal-approximation p-value.
    """
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return {"t": float("nan"), "df": float("nan"), "p_value": float("nan")}
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    denom = np.sqrt(vx / nx + vy / ny)
    if denom == 0:
        return {"t": float("nan"), "df": float("nan"), "p_value": 1.0}
    t = (mx - my) / denom
    df_num = (vx / nx + vy / ny) ** 2
    df_den = (vx ** 2) / (nx ** 2 * (nx - 1)) + (vy ** 2) / (ny ** 2 * (ny - 1))
    df = df_num / df_den if df_den > 0 else float("nan")
    try:
        from scipy.stats import t as t_dist
        p = 2.0 * (1.0 - t_dist.cdf(abs(t), df)) if np.isfinite(df) else float("nan")
    except ImportError:
        from math import erf, sqrt
        p = 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(t) / sqrt(2.0))))
    return {"t": float(t), "df": float(df), "p_value": float(p)}


def _mann_whitney(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """Two-sided Mann-Whitney U; returns NaN if scipy is unavailable."""
    try:
        from scipy.stats import mannwhitneyu
    except ImportError:
        return {"U": float("nan"), "p_value": float("nan")}
    if len(x) < 1 or len(y) < 1:
        return {"U": float("nan"), "p_value": float("nan")}
    try:
        res = mannwhitneyu(x, y, alternative="two-sided")
        return {"U": float(res.statistic), "p_value": float(res.pvalue)}
    except ValueError:
        return {"U": float("nan"), "p_value": float("nan")}


def _permutation_test(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 10_000,
    seed: Optional[int] = 0,
) -> Dict[str, float]:
    """
    Label-permutation two-sided test on the difference of means.

    Uses the exact enumeration when the number of partitions
    ``C(n_total, n_x)`` is <= ``n_perm``; otherwise falls back to a
    Monte-Carlo approximation with ``n_perm`` random relabellings.
    Returns the p-value and a flag indicating whether it was exact.
    """
    from math import comb
    nx, ny = len(x), len(y)
    if nx < 1 or ny < 1:
        return {"p_value": float("nan"), "n_perm": 0, "exact": False}
    combined = np.concatenate([x, y])
    n = nx + ny
    obs = abs(np.mean(x) - np.mean(y))

    total_combos = comb(n, nx)
    if total_combos <= n_perm:
        count = 0
        for idx in itertools.combinations(range(n), nx):
            mx = combined[list(idx)].mean()
            other = np.setdiff1d(np.arange(n), idx, assume_unique=True)
            my = combined[other].mean()
            if abs(mx - my) >= obs - 1e-12:
                count += 1
        return {"p_value": count / total_combos,
                "n_perm": int(total_combos), "exact": True}

    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        mx = combined[perm[:nx]].mean()
        my = combined[perm[nx:]].mean()
        if abs(mx - my) >= obs - 1e-12:
            count += 1
    return {"p_value": count / n_perm, "n_perm": int(n_perm), "exact": False}


def _cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    """Cohen's d with the pooled within-group std (ddof=1)."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan")
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if pooled == 0:
        return float("nan")
    return float((float(np.mean(x)) - float(np.mean(y))) / pooled)


def _summary(x: np.ndarray) -> Dict[str, float]:
    """Return (n, mean, std, sem) for a 1D array; NaNs produce None-like."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n == 0:
        return {"n": 0, "mean": float("nan"),
                "std": float("nan"), "sem": float("nan")}
    m = float(np.mean(x))
    if n == 1:
        return {"n": 1, "mean": m, "std": 0.0, "sem": 0.0}
    s = float(np.std(x, ddof=1))
    se = float(s / np.sqrt(n))
    return {"n": n, "mean": m, "std": s, "sem": se}


def _compare_vectors(
    x: np.ndarray,
    y: np.ndarray,
    name_x: str,
    name_y: str,
    log_transform: bool = False,
    n_boot: int = 10_000,
    seed: Optional[int] = 0,
) -> Dict[str, object]:
    """Pairwise comparison of two per-replicate vectors."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if log_transform:
        # Keep only strictly positive values for the log-space comparison.
        x = np.log10(x[x > 0])
        y = np.log10(y[y > 0])

    diff, ci_lo, ci_hi = _bootstrap_diff_ci(x, y, n_boot=n_boot, seed=seed)
    return {
        "condition_x": name_x,
        "condition_y": name_y,
        "log10_transformed": bool(log_transform),
        "summary_x": _summary(x),
        "summary_y": _summary(y),
        "diff_mean": diff,
        "diff_ci95": {"lo": ci_lo, "hi": ci_hi, "method": "percentile bootstrap",
                       "n_boot": int(n_boot)},
        "welch_t": _welch_t(x, y),
        "mann_whitney": _mann_whitney(x, y),
        "permutation": _permutation_test(x, y, n_perm=n_boot, seed=seed),
        "cohens_d": _cohens_d(x, y),
    }


def compare_vectors(
    x: Sequence[float],
    y: Sequence[float],
    name_x: str,
    name_y: str,
    log_transform: bool = False,
    n_boot: int = 10_000,
    seed: Optional[int] = 0,
) -> Dict[str, object]:
    """
    Public wrapper around ``_compare_vectors`` for reuse from other modules
    (e.g. ``scripts/summarize_analysis.py``). Runs the full between-
    condition battery: bootstrap CI on the mean difference, Welch's t,
    Mann-Whitney U, permutation, and Cohen's d, with optional log10
    transform for positive-valued, scale-spanning quantities.
    """
    return _compare_vectors(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        name_x, name_y,
        log_transform=log_transform,
        n_boot=n_boot, seed=seed,
    )


def compare_alpha(alpha_x: Sequence[float], alpha_y: Sequence[float],
                  name_x: str, name_y: str,
                  n_boot: int = 10_000, seed: Optional[int] = 0,
                  ) -> Dict[str, object]:
    """Wrapper for ``_compare_vectors`` on linear alpha values."""
    return _compare_vectors(np.asarray(alpha_x), np.asarray(alpha_y),
                             name_x, name_y, log_transform=False,
                             n_boot=n_boot, seed=seed)


def compare_K_alpha(K_x: Sequence[float], K_y: Sequence[float],
                    name_x: str, name_y: str,
                    n_boot: int = 10_000, seed: Optional[int] = 0,
                    ) -> Dict[str, object]:
    """
    Wrapper for ``_compare_vectors`` on K_alpha values.

    K_alpha is positive and varies over orders of magnitude, so the
    comparison is run in log10 space (as is standard for diffusion
    coefficients).  Effect sizes and CIs are reported in log10 units.
    """
    return _compare_vectors(np.asarray(K_x), np.asarray(K_y),
                             name_x, name_y, log_transform=True,
                             n_boot=n_boot, seed=seed)


# =============================================================================
# AGGREGATION FROM PER-REPLICATE JSON FILES
# =============================================================================

# Filename pattern produced by the orchestrator on copy-to-common:
#     <condition>_<nblocks>blk_rep<N>_msd_<label>.json
_PER_REP_MSD_RE = re.compile(
    r"^(?P<condition>.+)_(?P<nblocks>\d+)blk_rep(?P<rep>\d+)_msd_(?P<label>.+)\.json$"
)


def collect_per_replicate(common_dir: str) -> Dict[Tuple[str, str], List[dict]]:
    """
    Scan ``common_dir`` for per-replicate MSD JSONs and bucket them by
    ``(condition, pair_label)``.  Returns a dict whose values are lists
    of the loaded JSON payloads, one per replicate.

    A malformed or unreadable JSON is skipped with a warning.
    """
    buckets: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    pattern = os.path.join(common_dir, "*_msd_*.json")
    for path in sorted(glob.glob(pattern)):
        name = os.path.basename(path)
        m = _PER_REP_MSD_RE.match(name)
        if not m:
            continue
        key = (m.group("condition"), m.group("label"))
        try:
            with open(path) as f:
                payload = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"could not read {name}: {e}")
            continue
        payload["_source_file"] = name
        payload["_replicate"] = int(m.group("rep"))
        buckets[key].append(payload)
    return buckets


def _extract_per_rep_scalars(
    replicates: List[dict],
) -> Dict[str, np.ndarray]:
    """
    From a list of per-replicate MSD JSONs for the same (condition, pair),
    extract two vectors: alpha_per_rep and K_alpha_per_rep.  Missing fits
    are dropped.
    """
    alphas: List[float] = []
    Ks: List[float] = []
    for r in replicates:
        fit = r.get("alpha_fit") or {}
        a = fit.get("alpha")
        K = fit.get("K_alpha", fit.get("D"))
        if a is None or K is None:
            continue
        alphas.append(float(a))
        Ks.append(float(K))
    return {"alpha": np.array(alphas, dtype=float),
            "K_alpha": np.array(Ks, dtype=float)}


# =============================================================================
# ORCHESTRATOR
# =============================================================================

def run_msd_statistics(
    common_dir: str,
    out_json: str = "msd_between_condition_stats.json",
    out_forest_alpha: str = "msd_alpha_forest.png",
    out_forest_K: str = "msd_K_alpha_forest.png",
    n_boot: int = 10_000,
    seed: Optional[int] = 0,
    min_n_per_condition: int = 2,
) -> Optional[dict]:
    """
    Walk ``common_dir`` for per-replicate MSD JSONs, compute per-condition
    summaries (alpha and K_alpha), run pairwise comparisons between every
    pair of conditions that share the same probe pair, and emit::

        <common_dir>/<out_json>        : full summary payload
        <common_dir>/<out_forest_alpha>: forest plot of alpha across conditions
        <common_dir>/<out_forest_K>   : forest plot of log10(K_alpha)

    ``min_n_per_condition`` is the minimum number of replicates a
    condition must supply to be included in the comparisons.  The default
    of 2 keeps anything where at least one replicate fit succeeded on
    both sides of a comparison; raise to 3 in production.

    Returns the summary dict, or ``None`` if no replicates were found.
    """
    buckets = collect_per_replicate(common_dir)
    if not buckets:
        logger.info(f"  [msd-stats] no per-replicate MSD JSONs in {common_dir}")
        return None

    # Index by pair_label so we only compare conditions that used the same
    # probe pair.  With the conservation-aware auto-selection both cell
    # types should normally agree on the same 2 labels.
    by_pair: Dict[str, Dict[str, Dict[str, np.ndarray]]] = defaultdict(dict)
    for (condition, pair_label), reps in buckets.items():
        scalars = _extract_per_rep_scalars(reps)
        if scalars["alpha"].size < min_n_per_condition:
            logger.info(
                f"  [msd-stats] skipping ({condition}, {pair_label}): "
                f"only {scalars['alpha'].size} replicate fits")
            continue
        by_pair[pair_label][condition] = scalars

    summary: dict = {
        "schema_version": 1,
        "caveats": [
            "Statistics operate on per-replicate fits only.  Tiles within "
            "one replicate share a polymer chain and are not independent; "
            "they are reported as within-replicate uncertainty only.",
            "With n_rep = 3 per condition, Welch's t and permutation tests "
            "have limited power.  The smallest two-sided exact-permutation "
            "p-value achievable with 3 vs 3 is 2/20 = 0.10.",
            "Agreement across Welch, Mann-Whitney, and permutation is more "
            "informative than any single p-value.  Bootstrap CIs on the "
            "difference of means, and Cohen's d, quantify *how much* the "
            "conditions differ and should anchor interpretation.",
            "K_alpha is compared in log10 units (diffusion coefficients "
            "vary over orders of magnitude).",
        ],
        "n_boot_bootstrap": int(n_boot),
        "n_boot_permutation": int(n_boot),
        "pairs": {},
    }

    for pair_label, cond_map in sorted(by_pair.items()):
        pair_entry: dict = {
            "conditions": {},
            "comparisons_alpha": [],
            "comparisons_K_alpha": [],
        }
        for cond, sc in sorted(cond_map.items()):
            pair_entry["conditions"][cond] = {
                "alpha_per_rep": [float(x) for x in sc["alpha"]],
                "K_alpha_per_rep": [float(x) for x in sc["K_alpha"]],
                "alpha_summary": _summary(sc["alpha"]),
                "K_alpha_summary": _summary(sc["K_alpha"]),
                "log10_K_alpha_summary": _summary(
                    np.log10(sc["K_alpha"][sc["K_alpha"] > 0])),
            }

        conds_sorted = sorted(cond_map)
        for c1, c2 in itertools.combinations(conds_sorted, 2):
            a1, a2 = cond_map[c1]["alpha"], cond_map[c2]["alpha"]
            K1, K2 = cond_map[c1]["K_alpha"], cond_map[c2]["K_alpha"]
            pair_entry["comparisons_alpha"].append(
                compare_alpha(a1, a2, c1, c2, n_boot=n_boot, seed=seed)
            )
            pair_entry["comparisons_K_alpha"].append(
                compare_K_alpha(K1, K2, c1, c2, n_boot=n_boot, seed=seed)
            )
        summary["pairs"][pair_label] = pair_entry

    out_json_path = os.path.join(common_dir, out_json)
    with open(out_json_path, "w") as f:
        json.dump(summary, f, indent=2)
    logger.info(f"  [msd-stats] wrote {out_json_path}")

    # Forest plots: one panel per pair, one row per condition.
    try:
        for pair_label, pair_entry in summary["pairs"].items():
            plot_forest(
                pair_entry["conditions"],
                param="alpha",
                out_path=os.path.join(
                    common_dir, f"{pair_label}_{out_forest_alpha}"),
                title=f"alpha ({pair_label})",
                xlabel="anomalous exponent alpha",
            )
            plot_forest(
                pair_entry["conditions"],
                param="log10_K_alpha",
                out_path=os.path.join(
                    common_dir, f"{pair_label}_{out_forest_K}"),
                title=f"log10(K_alpha) ({pair_label})",
                xlabel="log10 K_alpha (monomer^2 / frame^alpha)",
            )
    except Exception as e:  # noqa: BLE001
        logger.warning(f"  [msd-stats] forest plot failed: {e}")

    return summary


# =============================================================================
# PLOTTING
# =============================================================================

def plot_forest(
    condition_entries: Dict[str, dict],
    param: str,
    out_path: str,
    title: str,
    xlabel: str,
    palette: str = "okabe-ito",
) -> None:
    """
    Forest plot of the per-replicate mean (+/- SE) of one parameter
    across conditions.  Dots = replicate values; solid bar = mean + SE;
    one row per condition.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if palette == "okabe-ito":
        colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7", "#000000"]
    else:
        cmap = plt.get_cmap(palette)
        colors = [cmap(i / max(1, len(condition_entries) - 1))
                  for i in range(len(condition_entries))]

    # Pick the per-replicate values to plot.
    values: Dict[str, List[float]] = {}
    for cond, entry in condition_entries.items():
        if param == "alpha":
            values[cond] = list(entry.get("alpha_per_rep", []))
        elif param == "log10_K_alpha":
            Ks = entry.get("K_alpha_per_rep", [])
            values[cond] = [float(np.log10(v)) for v in Ks if v and v > 0]
        else:
            raise ValueError(f"unknown param {param!r}")

    conds = list(values.keys())
    if not conds:
        return

    fig, ax = plt.subplots(figsize=(6.5, 0.5 + 0.45 * len(conds)))
    for i, cond in enumerate(conds):
        vals = np.asarray(values[cond], dtype=float)
        vals = vals[np.isfinite(vals)]
        y = i
        color = colors[i % len(colors)]
        if vals.size:
            m = float(np.mean(vals))
            se = float(np.std(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0
            ax.errorbar([m], [y], xerr=[[se], [se]], fmt="o",
                        color=color, ecolor=color, capsize=3, lw=1.6,
                        markersize=7)
            ax.plot(vals, np.full_like(vals, y, dtype=float), "o",
                    color=color, alpha=0.35, markersize=5)
    ax.set_yticks(range(len(conds)))
    ax.set_yticklabels(conds, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(True, axis="x", ls=":", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"    forest plot saved: {out_path}")
