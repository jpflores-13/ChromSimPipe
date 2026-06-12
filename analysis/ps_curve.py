#!/usr/bin/env python
"""
P(s) contact probability curve analysis — Hansen-lab-style.

This module is a thin, self-contained re-export of the P(s) machinery that
used to live inside contact_maps.py, plus a few convenience helpers for
per-condition overlay plotting and optional comparison with experimental
P(s) curves computed from a cooler file.

References
----------
Yang JH, Brandão HB, Hansen AS (2023) "DNA double-strand break end
synapsis by DNA loop extrusion", Nat. Commun. 14:1913,
DOI 10.1038/s41467-023-37583-w. Code:
  https://zenodo.org/records/7677969
  https://github.com/ahansenlab/DNA_break_synapsis_models
This is one of the Hansen-lab papers whose framing (loop extrusion +
log-log P(s) fitting) inspired our pipeline, but the P(s) helpers in
this module are NOT a direct port of any upstream file.

Hansen-lab P(s) reference implementation (renamed 2025–2026):
  https://github.com/ahansenlab/AbsQuant_analysis_code  (was AbsLoopQuant_analysis_code)
Their ``looptools.calculate_and_save_avg_Ps_curve`` wraps
``cooltools.expected_cis(smooth=True, aggregate_smoothed=True)`` over
all chromosomes; they extract no power-law exponent. The piecewise
log-log fit below (6th-degree polynomial below ``lower_bound_bp``,
1st-degree linear above) is the Mirny-lab / Gassler-style convention
we apply on the simulated dense maps; it is OUR design, not theirs.

Public API
----------
compute_ps_curve             : diagonal-average P(s) from a dense contact matrix
fit_ps_powerlaw              : piecewise polynomial + power-law log-log fit
extract_ps_metrics           : standardised summary metrics (P(s=X), alpha, slopes)
compute_ps_from_cooler       : Hansen-style cooltools wrapper (optional)
plot_ps_overlay              : overlay simulated (and optional experimental) P(s)
save_ps_json                 : persist metrics + fit to JSON

compute_ps_derivative        : d log P / d log s on a common log grid
plot_ps_derivative_overlay   : compare derivatives of several conditions
summarize_ps_derivatives     : one-row-per-condition table of derivative metrics
save_ps_derivative_table     : write that table as CSV + JSON
"""

from __future__ import annotations

import json
import logging
import os
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# CORE P(s) CALCULATION
# =============================================================================

def compute_ps_curve(
    contact_map: np.ndarray,
    resolution: int = 1000,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute P(s) — contact probability versus genomic distance — from a dense
    contact matrix by averaging each diagonal.

    Parameters
    ----------
    contact_map : np.ndarray
        Square contact frequency matrix.
    resolution : int
        Bp per bin.

    Returns
    -------
    distances : np.ndarray  (bp)
    ps        : np.ndarray  (contact probability at that distance)
    """
    N = contact_map.shape[0]
    ps = np.array([np.mean(np.diagonal(contact_map, offset=d)) for d in range(N)])
    distances = np.arange(N) * resolution
    return distances, ps


# =============================================================================
# POWER-LAW FITTING (Hansen lab style)
# =============================================================================

def fit_ps_powerlaw(
    distances: np.ndarray,
    ps: np.ndarray,
    lower_bound_bp: int = 2_000,
    upper_bound_bp: int = 30_000,
) -> Dict[str, object]:
    """
    Piecewise fit of P(s) in log-log space:
      - s < lower_bound_bp : 6th-degree polynomial (short-range, non-scaling)
      - lower_bound_bp <= s < upper_bound_bp : 1st-degree linear
        → exponent alpha (P(s) ~ s^-alpha).
    """
    valid = (distances > 0) & (ps > 0) & np.isfinite(ps)
    log_d = np.log10(distances[valid])
    log_p = np.log10(ps[valid])

    out: Dict[str, object] = {
        "lower_bound_bp": lower_bound_bp,
        "upper_bound_bp": upper_bound_bp,
    }

    d_valid = distances[valid]
    mask1 = d_valid < lower_bound_bp
    if mask1.sum() > 7:
        out["coeffs_region1"] = np.polyfit(log_d[mask1], log_p[mask1], deg=6).tolist()
    else:
        out["coeffs_region1"] = None

    mask2 = (d_valid >= lower_bound_bp) & (d_valid < upper_bound_bp)
    if mask2.sum() > 2:
        coeffs2 = np.polyfit(log_d[mask2], log_p[mask2], deg=1)
        out["coeffs_region2"] = coeffs2.tolist()
        out["exponent"] = float(-coeffs2[0])
        out["log_intercept"] = float(coeffs2[1])
    else:
        out["coeffs_region2"] = None
        out["exponent"] = None
        out["log_intercept"] = None

    return out


def extract_ps_metrics(
    distances: np.ndarray,
    ps: np.ndarray,
    resolution: int = 1000,
) -> Dict[str, float]:
    """
    Summary dict: P(s=10kb), ..., P(s=1Mb), power-law exponent, running
    slope at 5/10/50/100/500 kb. Missing entries are skipped silently.
    """
    metrics: Dict[str, float] = {}

    for label, target_bp in [
        ("10kb", 10_000), ("50kb", 50_000),
        ("100kb", 100_000), ("200kb", 200_000),
        ("500kb", 500_000), ("1Mb", 1_000_000),
    ]:
        idx = int(round(target_bp / resolution))
        if 0 <= idx < len(ps):
            metrics[f"P(s={label})"] = float(ps[idx])

    fit = fit_ps_powerlaw(distances, ps)
    if fit["exponent"] is not None:
        metrics["ps_exponent_2_30kb"] = float(fit["exponent"])
        metrics["ps_log_intercept"] = float(fit["log_intercept"])

    valid = (distances > 0) & (ps > 0) & np.isfinite(ps)
    d_v = distances[valid]
    p_v = ps[valid]
    if len(d_v) > 5:
        log_d = np.log10(d_v)
        log_p = np.log10(p_v)
        dlogp_dlogs = np.gradient(log_p, log_d)
        for label, target_bp in [
            ("5kb", 5_000), ("10kb", 10_000),
            ("50kb", 50_000), ("100kb", 100_000),
            ("500kb", 500_000),
        ]:
            idx = int(np.argmin(np.abs(d_v - target_bp)))
            if 2 <= idx < len(dlogp_dlogs) - 2:
                metrics[f"dlogP_dlogs_at_{label}"] = float(dlogp_dlogs[idx])

    return metrics


# =============================================================================
# EXPERIMENTAL P(s) FROM COOLER
# =============================================================================

def compute_ps_from_cooler(
    cool_path: str,
    region: Optional[str] = None,
    balance: bool = True,
    min_diag: int = 2,
    resolution: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute P(s) from a cooler file, restricted to a single cis region.

    If ``region`` is None the entire cis genome is used and the returned
    distances span the longest chromosome. Requires cooler + cooltools; if
    either is missing, raises ImportError so the caller can decide what to
    do (typically: skip the experimental comparison).

    Parameters
    ----------
    cool_path : str
        Path to an .mcool::/resolutions/<bp> URI or a single-resolution .cool.
    region : str, optional
        UCSC-style region, e.g. "chr3:34,000,000-36,000,000".
    balance : bool
        Use ICE-balanced values if available.
    min_diag : int
        Ignore the first ``min_diag`` diagonals (near-diagonal noise).
    resolution : int, optional
        Required if ``cool_path`` is an .mcool without an explicit resolution
        in the URI.

    Returns
    -------
    distances : np.ndarray (bp)
    ps        : np.ndarray
    """
    try:
        import cooler
    except ImportError as e:
        raise ImportError("cooler is required for compute_ps_from_cooler") from e

    # Resolve .mcool → single-resolution cooler
    uri = cool_path
    if cool_path.endswith(".mcool") and "::" not in cool_path:
        if resolution is None:
            raise ValueError(
                "resolution= must be provided for .mcool files without "
                "an explicit ::/resolutions/<bp> suffix"
            )
        uri = f"{cool_path}::/resolutions/{int(resolution)}"

    clr = cooler.Cooler(uri)
    bp_per_bin = clr.binsize

    if region is not None:
        mat = clr.matrix(balance=balance).fetch(region)
    else:
        # Largest cis block — pick the longest chromosome
        chromsizes = clr.chromsizes
        chrom = str(chromsizes.idxmax())
        mat = clr.matrix(balance=balance).fetch(chrom)

    mat = np.asarray(mat, dtype=float)
    mat[~np.isfinite(mat)] = 0.0

    N = mat.shape[0]
    ps_vals = []
    for d in range(N):
        diag = np.diagonal(mat, offset=d)
        finite = diag[np.isfinite(diag)]
        ps_vals.append(float(np.mean(finite)) if finite.size else 0.0)
    ps = np.asarray(ps_vals, dtype=float)

    if min_diag > 0:
        ps[:min_diag] = np.nan

    distances = np.arange(N) * bp_per_bin
    return distances, ps


# =============================================================================
# OVERLAY PLOTTING
# =============================================================================

def plot_ps_overlay(
    curves: Dict[str, Tuple[np.ndarray, np.ndarray]],
    out_path: str,
    title: str = "P(s) contact probability",
    lower_bound_bp: int = 2_000,
    upper_bound_bp: int = 30_000,
    xlim_bp: Tuple[int, int] = (1_000, 2_000_000),
    show_fit: bool = True,
) -> None:
    """
    Plot multiple P(s) curves on one log-log axis.

    Parameters
    ----------
    curves : dict
        ``{label: (distances_bp, ps)}``. Labels ending in "(exp)" are drawn
        with dashed lines to distinguish simulation from experiment.
    out_path : str
        PNG path.
    show_fit : bool
        Overlay the power-law fit from fit_ps_powerlaw for each simulated
        curve and annotate exponents.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    annot_lines: List[str] = []

    for label, (d, p) in curves.items():
        d = np.asarray(d)
        p = np.asarray(p)
        valid = (d > 0) & np.isfinite(p) & (p > 0)
        is_exp = label.endswith("(exp)")
        ls = "--" if is_exp else "-"
        line = ax.loglog(d[valid], p[valid], ls, label=label, lw=1.8)[0]

        if show_fit and not is_exp and valid.sum() > 5:
            fit = fit_ps_powerlaw(d, p, lower_bound_bp, upper_bound_bp)
            if fit["exponent"] is not None:
                x_fit = np.logspace(np.log10(lower_bound_bp),
                                    np.log10(upper_bound_bp), 50)
                slope, intercept = fit["coeffs_region2"]
                y_fit = 10 ** (slope * np.log10(x_fit) + intercept)
                ax.loglog(x_fit, y_fit, ":", color=line.get_color(), lw=1.0)
                annot_lines.append(f"{label}: α = {fit['exponent']:.2f}")

    ax.set_xlim(*xlim_bp)
    ax.set_xlabel("Genomic separation s (bp)")
    ax.set_ylabel("Contact probability P(s)")
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="lower left", fontsize=8)
    if annot_lines:
        ax.text(
            0.98, 0.98, "\n".join(annot_lines),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9),
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  P(s) overlay saved → {out_path}")


# =============================================================================
# PERSISTENCE HELPERS
# =============================================================================

def save_ps_json(
    out_path: str,
    distances: np.ndarray,
    ps: np.ndarray,
    resolution: int = 1000,
    extra: Optional[dict] = None,
) -> Dict:
    """Save P(s) summary metrics + fit to a JSON file and return the dict."""
    metrics = extract_ps_metrics(distances, ps, resolution)
    fit = fit_ps_powerlaw(distances, ps)
    payload = {"metrics": metrics, "fit": fit}
    if extra:
        payload.update(extra)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    return payload


def collect_curves(
    per_condition: Iterable[Tuple[str, np.ndarray]],
    resolution: int = 1000,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Convenience: ``[(label, contact_map), ...]`` → ``{label: (d, ps)}``.
    """
    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for label, cm in per_condition:
        d, p = compute_ps_curve(cm, resolution)
        out[label] = (d, p)
    return out


# =============================================================================
# DERIVATIVE OF P(s) ON LOG-LOG AXES
# =============================================================================

def _moving_average(y: np.ndarray, window: int) -> np.ndarray:
    """
    Centered moving average with edge reflection. Kept deliberately simple
    (no scipy dependency) so the module stays self-contained.
    """
    window = int(max(1, window))
    if window == 1 or y.size < 2 * window + 1:
        return y.copy()
    pad = window // 2
    pad_left = y[1:pad + 1][::-1]
    pad_right = y[-pad - 1:-1][::-1]
    padded = np.concatenate([pad_left, y, pad_right])
    kernel = np.ones(window) / window
    return np.convolve(padded, kernel, mode="valid")


def _gaussian_smooth_1d(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    """
    One-dimensional Gaussian smoothing in index space, matching the
    Gassler et al. 2017 convention (their Methods: scipy.ndimage.
    filters.gaussian_smoothing1d with radius ~ 0.8 after resampling onto
    a uniform log-s grid). Implemented with a reflected FIR kernel so
    there is no scipy dependency.
    """
    sigma = float(max(1e-6, sigma_bins))
    # Truncate at +- 4 sigma
    half = max(1, int(np.ceil(4.0 * sigma)))
    x = np.arange(-half, half + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()
    if y.size < 2 * half + 1:
        return y.copy()
    pad_left = y[1:half + 1][::-1]
    pad_right = y[-half - 1:-1][::-1]
    padded = np.concatenate([pad_left, y, pad_right])
    return np.convolve(padded, kernel, mode="valid")


def compute_ps_derivative(
    distances: np.ndarray,
    ps: np.ndarray,
    s_min_bp: float = 5_000.0,
    s_max_bp: float = 2_000_000.0,
    n_grid: int = 80,
    smooth_window: int = 5,
    smoothing: str = "gaussian",
    gaussian_sigma_bins: float = 0.8,
) -> Dict[str, np.ndarray]:
    """
    Compute the log-log slope (first derivative) of P(s):

        d log10(P)
        ----------        as a function of s, on a uniform log-s grid.
        d log10(s)

    The raw P(s) is interpolated into log space, lightly smoothed with a
    centred moving average of width ``smooth_window``, and differentiated
    with ``numpy.gradient``. The result is returned on a uniform log-s
    grid (``n_grid`` points between ``s_min_bp`` and ``s_max_bp``) so that
    multiple conditions can be compared pointwise.

    The derivative is the *local* contact-scaling exponent and is the
    quantity most directly interpretable for a mixer of polymer regimes:
    a flat region at ``-1.0`` to ``-1.5`` is Rouse/equilibrium-globule-like;
    a shallow region is loop-enriched; a steep region is TAD-boundary-like.

    Parameters
    ----------
    smoothing : {"gaussian", "moving_average"}
        Default ``"gaussian"`` matches the Gassler et al. 2017 convention
        (Methods: 1-D Gaussian smoothing of log-P and log-s with radius
        ``0.8`` index units). ``"moving_average"`` uses a centred box
        filter of width ``smooth_window`` and is kept for reproducibility
        with earlier runs.
    gaussian_sigma_bins : float
        Sigma of the Gaussian smoother in log-s grid-index units. Only
        used when ``smoothing == "gaussian"``. Default ``0.8`` follows
        Gassler 2017.

    Returns
    -------
    dict
        ``{"s": s_grid_bp, "dlogP_dlogs": slope, "logP": logP_on_grid}``
        All arrays length ``n_grid``.
    """
    d = np.asarray(distances, dtype=float)
    p = np.asarray(ps, dtype=float)

    valid = (d > 0) & (p > 0) & np.isfinite(p) & np.isfinite(d)
    if valid.sum() < 10:
        return {
            "s": np.array([]),
            "dlogP_dlogs": np.array([]),
            "logP": np.array([]),
        }

    log_d = np.log10(d[valid])
    log_p = np.log10(p[valid])

    # Restrict to the requested window and to where data exist
    s_lo = max(s_min_bp, d[valid].min())
    s_hi = min(s_max_bp, d[valid].max())
    if not (s_hi > s_lo):
        return {
            "s": np.array([]),
            "dlogP_dlogs": np.array([]),
            "logP": np.array([]),
        }

    s_grid = np.logspace(np.log10(s_lo), np.log10(s_hi), int(n_grid))
    log_s_grid = np.log10(s_grid)

    # Interpolate log(P) onto the log-s grid
    logP_on_grid = np.interp(log_s_grid, log_d, log_p)

    # Light smoothing to suppress diagonal-noise before differentiation.
    # Gaussian (default) follows Gassler et al. 2017 Methods; box filter is
    # kept as an alternative.
    if smoothing == "gaussian":
        logP_smoothed = _gaussian_smooth_1d(logP_on_grid, gaussian_sigma_bins)
    else:
        logP_smoothed = _moving_average(logP_on_grid, smooth_window)

    # Log-log gradient
    dlogP_dlogs = np.gradient(logP_smoothed, log_s_grid)

    return {
        "s": s_grid,
        "dlogP_dlogs": dlogP_dlogs,
        "logP": logP_on_grid,
        "smoothing": smoothing,
        "gaussian_sigma_bins": float(gaussian_sigma_bins),
    }


def summarize_ps_derivatives(
    curves: Dict[str, Tuple[np.ndarray, np.ndarray]],
    reference_separations_bp: Optional[List[float]] = None,
    s_min_bp: float = 5_000.0,
    s_max_bp: float = 2_000_000.0,
    n_grid: int = 80,
    smooth_window: int = 5,
    tad_window_bp: Tuple[float, float] = (100_000.0, 1_000_000.0),
    shortrange_window_bp: Tuple[float, float] = (5_000.0, 30_000.0),
    loop_size_search_window_bp: Tuple[float, float] = (20_000.0, 2_000_000.0),
) -> Dict[str, Dict[str, float]]:
    """
    One-row-per-condition table of descriptors of the P(s) derivative.

    For every condition in ``curves`` (``{label: (distances, ps)}``) compute
    the log-log derivative on a common grid and summarise it with a handful
    of interpretable numbers:

      - ``slope_at_<X>kb``     : d log P / d log s at the target separation
                                 (via linear interpolation on the log-s grid).
      - ``slope_min``, ``s_at_slope_min_bp`` : minimum (most negative) slope
                                 and where it occurs (TAD-scale "cliff").
      - ``slope_max``, ``s_at_slope_max_bp`` : maximum (least negative)
                                 slope and where it occurs (shoulder /
                                 loop-plateau signature).
      - ``slope_mean_shortrange`` / ``slope_mean_tadrange`` : mean slope in
                                 the short-range (default 5–30 kb) and
                                 TAD-range (default 100 kb–1 Mb) windows.

    Returns
    -------
    dict
        Nested dict: ``{condition_label: {metric_name: float}}``. Entries
        that cannot be evaluated are skipped silently.
    """
    if reference_separations_bp is None:
        reference_separations_bp = [10_000.0, 50_000.0, 100_000.0,
                                    500_000.0, 1_000_000.0]

    table: Dict[str, Dict[str, float]] = {}

    for label, (d, p) in curves.items():
        deriv = compute_ps_derivative(
            np.asarray(d), np.asarray(p),
            s_min_bp=s_min_bp, s_max_bp=s_max_bp,
            n_grid=n_grid, smooth_window=smooth_window,
        )
        s_grid = deriv["s"]
        slope = deriv["dlogP_dlogs"]
        if s_grid.size == 0:
            continue

        log_s_grid = np.log10(s_grid)
        row: Dict[str, float] = {}

        for s_ref in reference_separations_bp:
            if s_ref < s_grid[0] or s_ref > s_grid[-1]:
                continue
            slope_at_s = float(np.interp(np.log10(s_ref), log_s_grid, slope))
            key = f"slope_at_{_human_bp(s_ref)}"
            row[key] = slope_at_s

        i_min = int(np.argmin(slope))
        row["slope_min"] = float(slope[i_min])
        row["s_at_slope_min_bp"] = float(s_grid[i_min])

        # Gassler et al. 2017 interpretation: the location of the maximum
        # of d log P / d log s (smallest-in-magnitude slope) in the
        # interior of the curve gives the average cohesin-extruded loop
        # size. We restrict the argmax to an interior physical window
        # so the short-s edge of the grid (where the smoother has not
        # yet converged) cannot spuriously win.
        loop_mask = ((s_grid >= loop_size_search_window_bp[0]) &
                     (s_grid <= loop_size_search_window_bp[1]))
        if loop_mask.any():
            idx_in_window = np.where(loop_mask)[0]
            i_loopmax = int(idx_in_window[np.argmax(slope[loop_mask])])
            row["slope_max"] = float(slope[i_loopmax])
            row["s_at_slope_max_bp"] = float(s_grid[i_loopmax])
            row["inferred_mean_loop_size_bp"] = float(s_grid[i_loopmax])
            row["loop_size_search_window_bp"] = list(loop_size_search_window_bp)  # type: ignore[assignment]
        else:
            i_max = int(np.argmax(slope))
            row["slope_max"] = float(slope[i_max])
            row["s_at_slope_max_bp"] = float(s_grid[i_max])
            row["inferred_mean_loop_size_bp"] = float(s_grid[i_max])

        # The depth of the *minimum* at higher s (a proxy for the
        # linear density of loop-extruding cohesin in the same paper).
        row["loop_density_proxy_slope_min_depth"] = float(-slope[i_min])

        sr_mask = (s_grid >= shortrange_window_bp[0]) & (s_grid <= shortrange_window_bp[1])
        if sr_mask.any():
            row["slope_mean_shortrange"] = float(np.mean(slope[sr_mask]))
            row["slope_shortrange_window_bp"] = list(shortrange_window_bp)  # type: ignore[assignment]

        tad_mask = (s_grid >= tad_window_bp[0]) & (s_grid <= tad_window_bp[1])
        if tad_mask.any():
            row["slope_mean_tadrange"] = float(np.mean(slope[tad_mask]))
            row["slope_tad_window_bp"] = list(tad_window_bp)  # type: ignore[assignment]

        table[label] = row

    return table


def _human_bp(s_bp: float) -> str:
    """Format a separation in bp as '10kb', '1Mb', etc."""
    if s_bp >= 1_000_000:
        val = s_bp / 1_000_000
        return f"{val:.0f}Mb" if abs(val - round(val)) < 1e-6 else f"{val:g}Mb"
    if s_bp >= 1_000:
        val = s_bp / 1_000
        return f"{val:.0f}kb" if abs(val - round(val)) < 1e-6 else f"{val:g}kb"
    return f"{s_bp:.0f}bp"


def plot_ps_derivative_overlay(
    curves: Dict[str, Tuple[np.ndarray, np.ndarray]],
    out_path: str,
    title: str = "d log P(s) / d log s  (local contact-scaling exponent)",
    s_min_bp: float = 5_000.0,
    s_max_bp: float = 2_000_000.0,
    n_grid: int = 80,
    smooth_window: int = 5,
    reference_separations_bp: Optional[List[float]] = None,
    guide_slopes: Tuple[float, ...] = (-1.0, -1.5),
    loop_size_search_window_bp: Tuple[float, float] = (20_000.0, 2_000_000.0),
) -> None:
    """
    Plot the log-log derivative of P(s) for several conditions on the same
    axes. Vertical dashed lines mark the reference separations used by the
    summary table; horizontal dotted lines mark canonical polymer-regime
    slopes (``guide_slopes``, default Rouse-like ``-1.5`` and
    equilibrium-globule ``-1.0``).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if reference_separations_bp is None:
        reference_separations_bp = [10_000.0, 50_000.0, 100_000.0,
                                    500_000.0, 1_000_000.0]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    annot_lines: List[str] = []
    for label, (d, p) in curves.items():
        deriv = compute_ps_derivative(
            np.asarray(d), np.asarray(p),
            s_min_bp=s_min_bp, s_max_bp=s_max_bp,
            n_grid=n_grid, smooth_window=smooth_window,
        )
        if deriv["s"].size == 0:
            continue
        is_exp = label.endswith("(exp)")
        ls = "--" if is_exp else "-"
        line = ax.semilogx(deriv["s"], deriv["dlogP_dlogs"], ls, lw=1.8,
                           label=label)[0]

        # Gassler 2017 annotation: mark the maximum (= inferred loop
        # size) inside a physically-sensible window so short-s edge
        # effects of the smoother cannot claim to be loops.
        s_arr = deriv["s"]
        slope = deriv["dlogP_dlogs"]
        mask = ((s_arr >= loop_size_search_window_bp[0]) &
                (s_arr <= loop_size_search_window_bp[1]))
        if mask.any():
            idx_pool = np.where(mask)[0]
            i_max = int(idx_pool[np.argmax(slope[mask])])
        else:
            i_max = int(np.argmax(slope))
        s_peak = float(s_arr[i_max])
        y_peak = float(slope[i_max])
        ax.axvline(s_peak, color=line.get_color(), lw=0.5, ls="-", alpha=0.4)
        ax.plot([s_peak], [y_peak], marker="v",
                color=line.get_color(), ms=6, mec="black", mew=0.5)
        annot_lines.append(f"{label}: loop ≈ {_human_bp(s_peak)}")

    for s_ref in reference_separations_bp:
        if s_min_bp <= s_ref <= s_max_bp:
            ax.axvline(s_ref, color="0.7", lw=0.6, ls=":")

    for g in guide_slopes:
        ax.axhline(g, color="0.5", lw=0.6, ls=":")
        ax.text(s_max_bp, g, f" slope = {g:g}",
                ha="left", va="center", fontsize=7, color="0.4")

    ax.set_xlabel("Genomic separation s (bp)")
    ax.set_ylabel(r"$d \log_{10} P(s) / d \log_{10} s$")
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(loc="best", fontsize=8)

    if annot_lines:
        ax.text(
            0.02, 0.02, "\n".join(annot_lines) +
            "\n(Gassler 2017: argmax ≈ avg. cohesin loop size)",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85),
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  P(s) derivative overlay saved → {out_path}")


def save_ps_derivative_table(
    table: Dict[str, Dict[str, float]],
    csv_path: str,
    json_path: Optional[str] = None,
    reference_separations_bp: Optional[List[float]] = None,
) -> None:
    """
    Persist the summary table to a CSV (one row per condition, columns
    sorted for consistency) and optionally a JSON with identical content.

    Columns used in the CSV header (in order):

        condition,
        slope_at_10kb, slope_at_50kb, slope_at_100kb, slope_at_500kb,
        slope_at_1Mb,
        slope_mean_shortrange, slope_mean_tadrange,
        slope_min, s_at_slope_min_bp,
        slope_max, s_at_slope_max_bp
    """
    if reference_separations_bp is None:
        reference_separations_bp = [10_000.0, 50_000.0, 100_000.0,
                                    500_000.0, 1_000_000.0]

    ref_cols = [f"slope_at_{_human_bp(s)}" for s in reference_separations_bp]
    other_cols = [
        "slope_mean_shortrange", "slope_mean_tadrange",
        "slope_min", "s_at_slope_min_bp",
        "slope_max", "s_at_slope_max_bp",
        "inferred_mean_loop_size_bp",
        "loop_density_proxy_slope_min_depth",
    ]
    header = ["condition"] + ref_cols + other_cols

    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    with open(csv_path, "w") as f:
        f.write(",".join(header) + "\n")
        for label, row in table.items():
            cells = [label]
            for k in ref_cols + other_cols:
                v = row.get(k, "")
                if isinstance(v, float):
                    cells.append(f"{v:.4f}")
                else:
                    cells.append(str(v))
            f.write(",".join(cells) + "\n")

    if json_path is not None:
        # Drop the bare list columns (windows) from the JSON payload so it
        # parses as strict key→number; keep them alongside as `*_window_bp`.
        json_table: Dict[str, Dict[str, object]] = {}
        for label, row in table.items():
            json_table[label] = {k: v for k, v in row.items()}
        with open(json_path, "w") as f:
            json.dump({
                "metric_description": (
                    "d log10(P(s)) / d log10(s) summary. "
                    "Negative values: canonical. More negative = steeper "
                    "fall-off at that s. See README section 6."
                ),
                "reference_separations_bp": reference_separations_bp,
                "table": json_table,
            }, f, indent=2)

    logger.info(f"  P(s) derivative table saved → {csv_path}")
