#!/usr/bin/env python
"""
Two-point mean-squared displacement (MSD) analysis.

What this module computes
-------------------------
Given a simulated polymer trajectory (a sequence of (N_beads, 3) coordinate
arrays saved at regular intervals), this module measures how the vector
connecting two chosen monomers grows with time. Intuitively: if you stick
two fluorescent tags on a genome separated by some genomic distance and
watch them in a live cell, the two-point MSD tells you how their physical
separation fluctuates, and those fluctuations reveal whether the pair is
free-diffusing (MSD grows forever), trapped in a loop (MSD plateaus), or
somewhere in between.

Formally, for a pair (a, b) of monomer indices, define the separation
vector at frame t:

    r(t) = pos_b(t) : pos_a(t)

The time-averaged two-point MSD at lag tau is:

    MSD_2pt(tau) = < |r(t+tau) : r(t)|^2 >_t

which equals the MSD of the relative coordinate. By construction the pair's
common center-of-mass drift cancels, so the signal is not contaminated by
bulk translation. A Rouse (unlooped) polymer follows MSD_2pt ~ tau^{0.5};
a pair trapped in a cohesin-stabilised loop shows a plateau at 2*J (with J
the steady-state variance of the loop's end-to-end fluctuation).

Why we care for the cohesin-simulation project
----------------------------------------------
Contact maps are static summaries: they tell us how often two regions are
close, but not for how long or how fast they transition between states.
Live-cell imaging + MSD analysis in the Gabriele 2022 Fbn2 paper, and in
the Hansen-lab chromatin_dynamics package, was the direct experimental
handle that pinned down cohesin processivity, density, and capture rate.
If our simulations are right, the MSD curves from the simulation should
look like the ones from single-molecule tracking.

Key references
--------------
  Gabriele M et al. (2022) Science 376:496-501, Fbn2 two-color tracking.
  Mazzocca M, Narducci DN, Grosse-Holz S, Matthias J, Hansen AS (2025)
  bioRxiv 2025.05.10.653248, "Chromatin Dynamics are Highly Subdiffusive
  Across Seven Orders of Magnitude", code at
  github.com/ahansenlab/chromatin_dynamics. That preprint is the
  reference convention we follow for alpha / K_alpha:
  MSD(tau) = Gamma * tau^alpha, and bayesmsd uses the log-space
  parametrisation ``(log(alpha*Gamma), alpha)``. In our naming the
  prefactor ``Gamma`` becomes ``K_alpha``, so MSD(tau) = K_alpha * tau^alpha
  and K_alpha equals the MSD at tau = 1 saved-frame unit.
  Yang JH, Brandão HB, Hansen AS (2023) Nat Commun 14:1913,
  "DNA double-strand break end synapsis by DNA loop extrusion",
  DOI 10.1038/s41467-023-37583-w; code at
  github.com/ahansenlab/DNA_break_synapsis_models. Inspires the polymer
  simulation patterns; not a direct source for the MSD code below.
  Fbn2_simulations_and_data_analysis-main (local zip), dx_sim / dt_sim
  calibration formulas.

Simulation data has no localization error and no motion blur, so the
plain log-log polynomial fit below is appropriate. A Bayesian fit with
motion-blur priors (as in bayesmsd) is only needed when the trajectory
includes instrument noise; we keep the door open to plugging bayesmsd
in via ``fit_lag_min`` / ``fit_lag_max_frac`` if we later simulate
tagged-locus imaging.

Public API
----------
compute_two_point_msd        : MSD(tau) from a list of (N,3) frames.
compute_two_point_msd_tiled  : same, but extracts one pair per locus tile
                               for N_tiles-fold more statistics.
fit_msd_alpha                : log-log power-law fit over a chosen window.
fit_msd_saturation           : semi-empirical plateau + timescale fit.
pair_distance_distribution   : steady-state histogram of |r(t)|.
distance_autocorrelation     : < r(t)·r(t+tau) > / < |r|^2 > (loop memory).
plot_msd_curves              : overlay of simulated (and optional expt) MSD.
save_msd_json                : persist summary metrics to JSON.
apply_calibration            : convert sim units to nm and seconds.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# CORE MSD COMPUTATION
# =============================================================================

def _pair_trajectory(frames: Sequence[np.ndarray], idx_a: int, idx_b: int) -> np.ndarray:
    """
    Build the (n_frames, 3) separation-vector trajectory r(t) = pos_b - pos_a.

    Streams the input once, retaining only the two indexed positions per
    frame, so a full (T, N_chrom, 3) array is never allocated (would be
    168 GB for T=100k, N_chrom=70k).
    """
    T = len(frames)
    traj = np.empty((T, 3), dtype=np.float64)
    for t, frame in enumerate(frames):
        traj[t] = frame[idx_b] - frame[idx_a]
    return traj


def _pair_trajectories_tiled(
    frames: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    tile_size: int,
    pad: int,
    n_tiles: int,
) -> np.ndarray:
    """Per-tile separation-vector trajectories in a single streaming pass.

    Returns shape (n_tiles, T, 3). Callers that previously looped over
    ``range(n_tiles)`` calling ``_pair_trajectory`` (which re-read the
    HDF5 stream every tile) should use this instead.
    """
    tile_offsets = np.arange(n_tiles, dtype=int) * tile_size + pad
    a_idx = tile_offsets + idx_a
    b_idx = tile_offsets + idx_b
    T = len(frames)
    trajs = np.empty((n_tiles, T, 3), dtype=np.float64)
    for t, frame in enumerate(frames):
        trajs[:, t, :] = frame[b_idx] - frame[a_idx]
    return trajs


def _msd_from_trajectory(traj: np.ndarray,
                          lag_max: int,
                          lag_min: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Time-averaged MSD of a single separation-vector trajectory.

    MSD(tau) = < |traj[t+tau] - traj[t]|^2 >_t

    Uses O(T * lag_max) direct averaging. For the small lag range we
    actually fit (at most ~T/3 lags, trajectories of a few thousand frames),
    this is faster and much easier to debug than the FFT convolution path.

    Parameters
    ----------
    traj : np.ndarray of shape (T, 3)
        Separation-vector trajectory.
    lag_max : int
        Largest lag to include. Capped at T-1.
    lag_min : int
        Smallest lag (default 1).

    Returns
    -------
    lags : np.ndarray (int, length L)
    msd  : np.ndarray (float, length L)
    """
    T = traj.shape[0]
    lag_max = min(int(lag_max), T - 1)
    lag_min = max(1, int(lag_min))
    if lag_max < lag_min:
        return np.empty(0, dtype=int), np.empty(0, dtype=float)
    lags = np.arange(lag_min, lag_max + 1, dtype=int)
    msd = np.empty(lags.shape, dtype=np.float64)
    for k, tau in enumerate(lags):
        diff = traj[tau:] - traj[:-tau]
        sq = np.einsum("ij,ij->i", diff, diff)   # faster than (**2).sum(1)
        msd[k] = sq.mean()
    return lags, msd


def compute_two_point_msd(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    lag_min: int = 1,
    lag_max: Optional[int] = None,
    lag_max_frac: float = 1 / 3,
) -> Dict[str, np.ndarray]:
    """
    Compute the two-point MSD for a pair (idx_a, idx_b) from a single
    trajectory. Returns a dict with ``lags`` (int frames), ``msd`` (in sim
    length^2 units), and ``traj`` (separation vector trajectory, for
    downstream fits or localisation-error models).

    Parameters
    ----------
    conformations : list of np.ndarray
        Each (N_beads, 3). Consecutive entries represent successive saved
        simulation frames.
    idx_a, idx_b : int
        Monomer indices to track.
    lag_min, lag_max : int
        Lag range in frames. If ``lag_max`` is None, defaults to
        ``int(n_frames * lag_max_frac)``.

    Returns
    -------
    dict with keys: lags, msd, traj
    """
    if len(conformations) < 2:
        raise ValueError("Need at least 2 conformations to compute an MSD.")

    traj = _pair_trajectory(conformations, idx_a, idx_b)

    if lag_max is None:
        lag_max = max(1, int(traj.shape[0] * lag_max_frac))

    lags, msd = _msd_from_trajectory(traj, lag_max=lag_max, lag_min=lag_min)
    return {"lags": lags, "msd": msd, "traj": traj}


def compute_two_point_msd_tiled(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    tile_size: int,
    pad: int,
    n_tiles: int,
    lag_min: int = 1,
    lag_max: Optional[int] = None,
    lag_max_frac: float = 1 / 3,
) -> Dict[str, np.ndarray]:
    """
    Tile-aware wrapper around :func:`compute_two_point_msd`.

    The simulation uses the Hansen-lab tiling trick: one long chromosome
    contains ``n_tiles`` identical copies of the 2 Mb locus, each offset by
    ``tile_size`` monomers and surrounded by a ``pad`` of padding. For a
    pair given in *locus* coordinates (0..N_locus-1), this function
    extracts the corresponding pair in every tile, computes the MSD per
    tile, and averages the MSDs (equivalent to concatenating independent
    two-point trajectories).

    Parameters
    ----------
    conformations : list of np.ndarray
        Each (chrom_size, 3), where chrom_size == n_tiles * tile_size.
    idx_a, idx_b : int
        Monomer indices in locus coordinates (0..tile_size - 2*pad - 1).
    tile_size, pad, n_tiles : int
        From ``configs.parameters.TILING``.

    Returns
    -------
    dict: lags, msd, n_tiles_used, per_tile_msd (n_tiles, L).
        ``msd`` is the mean across tiles; ``per_tile_msd`` keeps each tile's
        curve for uncertainty bands.
    """
    # Accept either a list/tuple of per-frame arrays or a single stacked
    # (n_frames, n_beads, 3) array.  Use len() so numpy arrays work too.
    if len(conformations) == 0:
        raise ValueError("Empty conformation list.")
    chrom_size = conformations[0].shape[0]

    # Non-tiled conformations: fall back to a single trajectory
    if chrom_size == tile_size - 2 * pad or chrom_size < n_tiles * tile_size:
        # Either single locus or arbitrary single chromosome; don't tile.
        if idx_a >= chrom_size or idx_b >= chrom_size:
            raise ValueError(
                f"Pair ({idx_a}, {idx_b}) out of bounds for chromosome "
                f"size {chrom_size}."
            )
        return {
            **compute_two_point_msd(
                conformations, idx_a, idx_b,
                lag_min=lag_min, lag_max=lag_max, lag_max_frac=lag_max_frac,
            ),
            "n_tiles_used": 1,
            "per_tile_msd": None,
        }

    # Single streaming pass across all tiles.
    trajs_stack = _pair_trajectories_tiled(
        conformations, idx_a, idx_b, tile_size, pad, n_tiles,
    )                                                       # (n_tiles, T, 3)
    per_tile_trajs: List[np.ndarray] = [trajs_stack[k] for k in range(n_tiles)]

    T = per_tile_trajs[0].shape[0]
    if lag_max is None:
        lag_max = max(1, int(T * lag_max_frac))

    per_tile_msd = []
    lags_ref: Optional[np.ndarray] = None
    for traj in per_tile_trajs:
        lags, msd = _msd_from_trajectory(traj, lag_max=lag_max, lag_min=lag_min)
        per_tile_msd.append(msd)
        lags_ref = lags if lags_ref is None else lags_ref
    per_tile_msd_arr = np.vstack(per_tile_msd)
    mean_msd = per_tile_msd_arr.mean(axis=0)

    return {
        "lags": lags_ref if lags_ref is not None else np.empty(0, dtype=int),
        "msd": mean_msd,
        "n_tiles_used": len(per_tile_trajs),
        "per_tile_msd": per_tile_msd_arr,
    }


# =============================================================================
# FITTING
# =============================================================================

def fit_msd_alpha(
    lags: np.ndarray,
    msd: np.ndarray,
    fit_lag_min: int = 5,
    fit_lag_max_frac: float = 0.25,
    min_n_lags_for_fit: int = 6,
) -> Dict[str, Optional[float]]:
    """
    Log-log linear fit to MSD in the sub-diffusive window.

    A Rouse polymer has MSD_2pt(tau) ~ tau^{0.5}; cohesin stabilisation
    bends the curve toward smaller alpha.  The functional form is::

        MSD(tau) = K_alpha * tau ** alpha

    so a straight line in log-log coordinates with slope ``alpha`` and
    intercept ``log10(K_alpha)`` recovers both parameters simultaneously.

    ``K_alpha`` is the generalised (anomalous) diffusion coefficient: it
    is the value of MSD at tau = 1 (in saved-frame units) and collapses to
    the classical diffusion coefficient only when alpha == 1.  Two
    conditions with the same alpha but different K_alpha differ in the
    *amplitude* of the pair fluctuation; two conditions with the same
    K_alpha but different alpha differ in the *memory* of the fluctuation
    (sub-diffusive vs. Rouse vs. super-diffusive).

    The returned dict keeps ``D`` as an alias of ``K_alpha`` for backward
    compatibility with older callers.

    Parameters
    ----------
    fit_lag_min : int
        Lower cutoff on lags (drop short-time bumps from bond relaxation).
    fit_lag_max_frac : float
        Upper cutoff as a fraction of max lag (drop long-time noise).
    min_n_lags_for_fit : int
        Minimum number of lags in the window to attempt the fit.

    Returns
    -------
    dict with keys:
        alpha             : float, slope of log-log fit (the anomalous exponent).
        K_alpha           : float, prefactor so MSD = K_alpha * tau^alpha.
        D                 : float, alias of K_alpha (legacy).
        log_intercept     : float, log10(K_alpha).
        alpha_stderr      : float, 1-sigma uncertainty on alpha from the
                            residual covariance of the log-log fit.
        log_K_alpha_stderr: float, 1-sigma uncertainty on log10(K_alpha).
        fit_lag_min, fit_lag_max, n_points.
        Missing entries are ``None`` if the window is too short.
    """
    lags = np.asarray(lags, dtype=float)
    msd = np.asarray(msd, dtype=float)
    empty = {"alpha": None, "K_alpha": None, "D": None,
             "log_intercept": None,
             "alpha_stderr": None, "log_K_alpha_stderr": None,
             "fit_lag_min": fit_lag_min, "fit_lag_max": None, "n_points": 0}
    if lags.size == 0:
        return empty

    fit_lag_max = max(fit_lag_min + 1, int(lags.max() * fit_lag_max_frac))
    mask = (
        (lags >= fit_lag_min) & (lags <= fit_lag_max)
        & (msd > 0) & np.isfinite(msd)
    )
    if mask.sum() < min_n_lags_for_fit:
        out = dict(empty)
        out["fit_lag_max"] = fit_lag_max
        out["n_points"] = int(mask.sum())
        return out

    log_t = np.log10(lags[mask])
    log_msd = np.log10(msd[mask])
    # polyfit with cov=True returns the covariance matrix of the fit
    # parameters, from which we can extract the 1-sigma SE on slope and
    # intercept.  This is purely the fit's internal uncertainty (how well
    # the points line up); it does NOT capture between-replicate variance.
    try:
        (slope, intercept), cov = np.polyfit(log_t, log_msd, deg=1, cov=True)
        alpha_se = float(np.sqrt(max(cov[0, 0], 0.0)))
        logk_se = float(np.sqrt(max(cov[1, 1], 0.0)))
    except (np.linalg.LinAlgError, ValueError):
        slope, intercept = np.polyfit(log_t, log_msd, deg=1)
        alpha_se = None
        logk_se = None
    K_alpha = float(10 ** intercept)
    return {
        "alpha": float(slope),
        "K_alpha": K_alpha,
        "D": K_alpha,                        # legacy alias: MSD = D * tau^alpha
        "log_intercept": float(intercept),
        "alpha_stderr": alpha_se,
        "log_K_alpha_stderr": logk_se,
        "fit_lag_min": int(fit_lag_min),
        "fit_lag_max": int(fit_lag_max),
        "n_points": int(mask.sum()),
    }


def fit_msd_alpha_per_tile(
    lags: np.ndarray,
    per_tile_msd: np.ndarray,
    fit_lag_min: int = 5,
    fit_lag_max_frac: float = 0.25,
    min_n_lags_for_fit: int = 6,
) -> Dict[str, Optional[np.ndarray]]:
    """
    Run :func:`fit_msd_alpha` on each tile's MSD curve independently.

    With the Hansen-lab tiling trick, one simulation yields ``n_tiles``
    independent copies of the two-point MSD for the same pair.  Fitting
    alpha and K_alpha on each tile produces an empirical distribution of
    the two parameters that captures the spread within one replicate.

    Caveat: tiles share a single polymer chain and a single LEF pool, so
    their per-tile fits are not strictly independent.  The spread is a
    reasonable *lower bound* on the uncertainty, but between-replicate
    variance is a stricter (and more appropriate) reference.  Use the
    returned arrays for within-replicate diagnostics, not for between-
    condition null-hypothesis tests.

    Parameters
    ----------
    lags : np.ndarray, shape (L,)
        Common lag axis (one vector shared by all tiles).
    per_tile_msd : np.ndarray, shape (n_tiles, L)
        One MSD curve per tile (from ``compute_two_point_msd_tiled``).
    fit_lag_min, fit_lag_max_frac, min_n_lags_for_fit
        Same as :func:`fit_msd_alpha`.

    Returns
    -------
    dict with:
        alpha_per_tile   : np.ndarray, shape (n_tiles,)
        K_alpha_per_tile : np.ndarray, shape (n_tiles,)
        alpha_mean, alpha_sem, alpha_std
        K_alpha_mean, K_alpha_sem, K_alpha_std
        n_tiles_fitted   : int (tiles for which the fit succeeded)
    """
    if per_tile_msd is None:
        return {"alpha_per_tile": None, "K_alpha_per_tile": None,
                "alpha_mean": None, "alpha_sem": None, "alpha_std": None,
                "K_alpha_mean": None, "K_alpha_sem": None, "K_alpha_std": None,
                "n_tiles_fitted": 0}
    per_tile_msd = np.asarray(per_tile_msd)
    if per_tile_msd.ndim != 2:
        raise ValueError(
            f"per_tile_msd must be (n_tiles, L), got shape {per_tile_msd.shape}")
    n_tiles = per_tile_msd.shape[0]
    alphas = np.full(n_tiles, np.nan, dtype=float)
    Ks = np.full(n_tiles, np.nan, dtype=float)
    for k in range(n_tiles):
        fit = fit_msd_alpha(
            lags, per_tile_msd[k],
            fit_lag_min=fit_lag_min,
            fit_lag_max_frac=fit_lag_max_frac,
            min_n_lags_for_fit=min_n_lags_for_fit,
        )
        if fit.get("alpha") is not None:
            alphas[k] = fit["alpha"]
            Ks[k] = fit["K_alpha"]
    valid = np.isfinite(alphas) & np.isfinite(Ks)
    n_valid = int(valid.sum())

    def _m_s_se(x):
        if n_valid == 0:
            return (None, None, None)
        xv = x[valid]
        m = float(np.mean(xv))
        s = float(np.std(xv, ddof=1)) if n_valid > 1 else 0.0
        se = float(s / np.sqrt(n_valid)) if n_valid > 1 else 0.0
        return (m, s, se)

    a_m, a_s, a_se = _m_s_se(alphas)
    K_m, K_s, K_se = _m_s_se(Ks)
    return {
        "alpha_per_tile": alphas,
        "K_alpha_per_tile": Ks,
        "alpha_mean": a_m, "alpha_std": a_s, "alpha_sem": a_se,
        "K_alpha_mean": K_m, "K_alpha_std": K_s, "K_alpha_sem": K_se,
        "n_tiles_fitted": n_valid,
    }


def fit_msd_saturation(
    lags: np.ndarray,
    msd: np.ndarray,
    loc_error_sq: float = 0.0,
) -> Dict[str, Optional[float]]:
    """
    Fit a semi-empirical saturating form::

        MSD(tau) = 2 * J * (1 - exp(-(tau / tau_c)^alpha)) + 2 * sigma^2

    This is a convenience shape: J is the "plateau half-height" (steady-
    state half-variance of the separation vector) and tau_c is the
    relaxation time to reach it. ``sigma^2`` is an optional localisation-
    error term (pass 0 for pure sim data). Robust to short or noisy
    trajectories by clipping MSD to 1e-12 before log transforms.

    Returns None-filled dict if the trajectory is too short or scipy
    is not available.
    """
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        return {"J": None, "tau_c": None, "alpha_sat": None, "sigma2": None}

    lags = np.asarray(lags, dtype=float)
    msd = np.asarray(msd, dtype=float)
    if lags.size < 8 or msd.size < 8:
        return {"J": None, "tau_c": None, "alpha_sat": None, "sigma2": None}

    def _shape(tau, J, tau_c, alpha):
        return 2.0 * J * (1.0 - np.exp(-((tau / max(tau_c, 1e-6)) ** alpha))) \
            + 2.0 * loc_error_sq

    try:
        p0 = [float(np.nanmax(msd)) / 2.0,
              float(lags[len(lags) // 3]),
              0.5]
        popt, _ = curve_fit(
            _shape, lags, msd, p0=p0,
            bounds=([0, 1e-3, 0.05], [np.inf, np.inf, 2.0]),
            maxfev=4000,
        )
        return {
            "J": float(popt[0]),
            "tau_c": float(popt[1]),
            "alpha_sat": float(popt[2]),
            "sigma2": float(loc_error_sq),
        }
    except Exception as e:  # noqa: BLE001
        logger.debug(f"fit_msd_saturation failed: {e}")
        return {"J": None, "tau_c": None, "alpha_sat": None, "sigma2": None}


# =============================================================================
# COMPLEMENTARY METRICS (pair distance distribution + autocorrelation)
# =============================================================================

def pair_distance_distribution(
    traj: np.ndarray,
    bins: int = 60,
    range_: Optional[Tuple[float, float]] = None,
) -> Dict[str, np.ndarray]:
    """
    Steady-state histogram of |r(t)|, the pair-separation magnitude.

    Useful to visualise whether the pair has a bimodal distribution
    (looped vs unlooped) or a broad unimodal one.
    """
    d = np.sqrt(np.einsum("ij,ij->i", traj, traj))
    counts, edges = np.histogram(d, bins=bins, range=range_, density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return {"centers": centers, "density": counts, "mean_d": float(d.mean()),
            "median_d": float(np.median(d)), "std_d": float(d.std())}


def distance_autocorrelation(
    traj: np.ndarray,
    lag_max: Optional[int] = None,
    lag_max_frac: float = 1 / 3,
) -> Dict[str, np.ndarray]:
    """
    Normalised autocorrelation of the separation-vector norm::

        C(tau) = < (|r(t)|-<|r|>) * (|r(t+tau)|-<|r|>) >_t / Var(|r|)

    A slower decay means the pair spends longer in the same conformation
    (e.g., sitting inside a persistent cohesin loop).
    """
    d = np.sqrt(np.einsum("ij,ij->i", traj, traj))
    d0 = d - d.mean()
    T = d0.size
    if lag_max is None:
        lag_max = max(1, int(T * lag_max_frac))
    var = d0.var()
    if var <= 0:
        return {"lags": np.arange(lag_max + 1, dtype=int),
                "C": np.zeros(lag_max + 1, dtype=float)}
    lags = np.arange(lag_max + 1, dtype=int)
    C = np.empty(lags.shape, dtype=float)
    for k, tau in enumerate(lags):
        if tau == 0:
            C[k] = 1.0
        else:
            C[k] = float(np.mean(d0[:-tau] * d0[tau:])) / var
    return {"lags": lags, "C": C}


# =============================================================================
# PHYSICAL-UNITS CALIBRATION
# =============================================================================

@dataclass
class Calibration:
    """
    Calibration from simulation units to (nm, seconds).

    Attributes
    ----------
    nm_per_monomer : float
        Effective physical size of one monomer in the coarse-grained polymer.
    sec_per_frame : float
        Real time elapsed between two saved simulation frames.
    source : str
        "hic" (Ps(s) match) or "msd" (microscopy MSD match) or "assumed".
    """
    nm_per_monomer: float
    sec_per_frame: float
    source: str

    def to_dict(self) -> dict:
        return asdict(self)


def calibrate_from_hic(
    bp_per_monomer: int,
    persistence_nm: float = 50.0,
    sec_per_frame: Optional[float] = None,
) -> Calibration:
    """
    Simplest calibration: assume 1 monomer = ``bp_per_monomer`` base pairs,
    and convert genomic length to physical length via a persistence-length
    estimate for chromatin (Gabriele 2022 used ~50 nm for 1 kb monomers of
    mESC chromatin).

    If ``sec_per_frame`` is not supplied, caller should plug one in later
    based on simulation block size and the extrusion rate (~1 kb/s).

    Returns
    -------
    Calibration
    """
    nm_per_mon = persistence_nm * np.sqrt(bp_per_monomer / 1000.0)
    return Calibration(
        nm_per_monomer=float(nm_per_mon),
        sec_per_frame=float(sec_per_frame) if sec_per_frame else 1.0,
        source="hic",
    )


def calibrate_from_msd(
    sim_lags_frames: np.ndarray,
    sim_msd_monomer2: np.ndarray,
    expt_lags_sec: np.ndarray,
    expt_msd_um2: np.ndarray,
    fit_lag_min_sec: float = 1.0,
    fit_lag_max_sec: float = 30.0,
) -> Calibration:
    """
    Match simulation MSD to experimental microscopy MSD.

    Implements the Fbn2 / calibrate_simulation_units_to_real_units.py
    approach::

        dx_sim = sqrt(J_expt / J_sim) * dx_expt
        dt_sim = ((G_sim * dx_sim^2) / (G_expt * dx_expt^2))^2 * dt_expt

    where J is the steady-state plateau and G is the short-time
    diffusivity. Here we take ``dx_expt = 1 um`` and ``dt_expt = 1 s`` so
    the returned units are um/monomer and s/frame; a post-hoc conversion
    to nm is applied inside the Calibration.

    Parameters
    ----------
    sim_lags_frames : np.ndarray
        Simulation lag axis, in frame units.
    sim_msd_monomer2 : np.ndarray
        Simulation MSD in monomer^2.
    expt_lags_sec, expt_msd_um2 : np.ndarray
        Experimental microscopy curve.
    fit_lag_min_sec, fit_lag_max_sec : float
        Window for G (diffusivity) fit on the experimental curve.
    """
    # Short-time diffusivities from log-log slope prefactors
    def _short_time_diffusivity(x, y, x_min, x_max):
        mask = (x >= x_min) & (x <= x_max) & (y > 0) & np.isfinite(y)
        if mask.sum() < 3:
            return np.nan
        slope, intercept = np.polyfit(np.log10(x[mask]), np.log10(y[mask]), deg=1)
        return float(10 ** intercept), float(slope)

    sim_valid = (sim_lags_frames > 0) & (sim_msd_monomer2 > 0)
    expt_valid = (expt_lags_sec > 0) & (expt_msd_um2 > 0)

    J_sim = float(np.nanmax(sim_msd_monomer2[sim_valid])) / 2.0
    J_expt = float(np.nanmax(expt_msd_um2[expt_valid])) / 2.0

    G_expt, _ = _short_time_diffusivity(expt_lags_sec, expt_msd_um2,
                                        fit_lag_min_sec, fit_lag_max_sec)
    # For sim, scale the experimental window to frames assuming
    # sec_per_frame=1 initially — then we iterate on the conversion.
    G_sim, _ = _short_time_diffusivity(sim_lags_frames, sim_msd_monomer2,
                                        fit_lag_min_sec, fit_lag_max_sec)

    dx_sim_um = float(np.sqrt(max(J_expt, 1e-12) / max(J_sim, 1e-12)))
    dt_sim_sec = float(((G_sim * dx_sim_um ** 2) /
                        max(G_expt * 1.0 ** 2, 1e-12)) ** 2)
    return Calibration(
        nm_per_monomer=dx_sim_um * 1000.0,
        sec_per_frame=dt_sim_sec,
        source="msd",
    )


def apply_calibration(
    lags_frames: np.ndarray,
    msd_monomer2: np.ndarray,
    calib: Calibration,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert MSD in (monomer, frame) units to (nm, second) units."""
    lags_sec = np.asarray(lags_frames, dtype=float) * calib.sec_per_frame
    msd_nm2 = np.asarray(msd_monomer2, dtype=float) * (calib.nm_per_monomer ** 2)
    return lags_sec, msd_nm2


# =============================================================================
# PLOTTING
# =============================================================================

def plot_msd_curves(
    curves: Dict[str, Dict[str, np.ndarray]],
    out_path: str,
    title: str = "Two-point MSD",
    xlabel: str = "lag τ (frames)",
    ylabel: str = "MSD (monomer²)",
    show_alpha: bool = True,
    rouse_guide: bool = True,
    fit_range: Optional[Tuple[int, int]] = None,
    palette: str = "okabe-ito",
) -> None:
    """
    Overlay multiple MSD curves (one per condition / pair) on a log-log axis.

    Parameters
    ----------
    curves : dict
        ``{label: {"lags": np.ndarray, "msd": np.ndarray,
                   "alpha": Optional[float], "fit_range": Optional[(lo, hi)]}}``.
    show_alpha : bool
        Print fitted alpha next to each label in the legend.
    rouse_guide : bool
        Overlay a dashed reference line MSD ~ tau^0.5 for comparison.
    palette : str
        ``"okabe-ito"`` (default, CVD-safe), ``"viridis"`` or ``"cividis"``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if palette == "okabe-ito":
        colors = [
            "#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#000000",
        ]
    else:
        cmap = plt.get_cmap(palette)
        colors = [cmap(i / max(1, len(curves) - 1)) for i in range(len(curves))]

    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    annot_lines: List[str] = []

    all_lags = []
    for i, (label, c) in enumerate(curves.items()):
        lags = np.asarray(c["lags"], dtype=float)
        msd = np.asarray(c["msd"], dtype=float)
        valid = (lags > 0) & (msd > 0) & np.isfinite(msd)
        color = colors[i % len(colors)]
        line = ax.loglog(lags[valid], msd[valid], "-", color=color,
                         label=label, lw=1.8)[0]
        all_lags.append(lags[valid])

        alpha = c.get("alpha")
        if show_alpha and alpha is not None:
            annot_lines.append(f"{label}: α = {alpha:.2f}")
            lo, hi = c.get("fit_range", (None, None))
            if lo is not None and hi is not None:
                x = np.logspace(np.log10(lo), np.log10(hi), 40)
                D = c.get("D", None)
                if D is not None:
                    y = D * x ** alpha
                    ax.loglog(x, y, ":", color=line.get_color(), lw=1.0)

    if rouse_guide and all_lags:
        lag_union = np.concatenate(all_lags)
        if lag_union.size:
            t = np.logspace(np.log10(max(1.0, lag_union.min())),
                            np.log10(lag_union.max()), 50)
            # Anchor the guide at the first finite sim MSD value to keep
            # it visually near (not on top of) the data.
            first_label = next(iter(curves))
            msd0 = np.asarray(curves[first_label]["msd"], dtype=float)
            lags0 = np.asarray(curves[first_label]["lags"], dtype=float)
            ok = (lags0 > 0) & (msd0 > 0) & np.isfinite(msd0)
            if ok.any():
                ref_lag = lags0[ok][0]
                ref_msd = msd0[ok][0]
                y = ref_msd * (t / ref_lag) ** 0.5
                ax.loglog(t, y, "--", color="0.4", lw=1.0,
                          label="Rouse: MSD ~ τ^0.5")

    if fit_range:
        ax.axvspan(fit_range[0], fit_range[1], color="0.9", alpha=0.3, zorder=0)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)

    if annot_lines:
        ax.text(0.98, 0.02, "\n".join(annot_lines),
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, family="monospace",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  MSD overlay saved: {out_path}")


def plot_pair_distance_distribution(
    per_label: Dict[str, Dict[str, np.ndarray]],
    out_path: str,
    title: str = "Pair distance distribution (steady state)",
    palette: str = "okabe-ito",
) -> None:
    """Overlay pair-separation histograms from multiple labels."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if palette == "okabe-ito":
        colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7", "#000000"]
    else:
        cmap = plt.get_cmap(palette)
        colors = [cmap(i / max(1, len(per_label) - 1))
                  for i in range(len(per_label))]

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for i, (label, d) in enumerate(per_label.items()):
        ax.plot(d["centers"], d["density"], "-", lw=1.6,
                color=colors[i % len(colors)], label=label)
    ax.set_xlabel("|r| (monomer units)")
    ax.set_ylabel("density")
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  Pair-distance density saved: {out_path}")


# =============================================================================
# PERSISTENCE
# =============================================================================

def save_msd_json(
    out_path: str,
    lags: np.ndarray,
    msd: np.ndarray,
    pair: Tuple[int, int],
    label: str,
    fit: Dict[str, Optional[float]],
    sat_fit: Optional[Dict[str, Optional[float]]] = None,
    distance_stats: Optional[Dict[str, float]] = None,
    calib: Optional[Calibration] = None,
    extra: Optional[dict] = None,
) -> Dict:
    """Write a compact JSON summary of one MSD measurement."""
    payload = {
        "label": label,
        "pair": {"monomer_a": int(pair[0]), "monomer_b": int(pair[1]),
                 "sep_monomers": abs(int(pair[1]) - int(pair[0]))},
        "n_lags": int(len(lags)),
        "lags_frames": [int(x) for x in lags],
        "msd_monomer2": [float(x) for x in msd],
        "alpha_fit": fit,
    }
    if sat_fit is not None:
        payload["saturation_fit"] = sat_fit
    if distance_stats is not None:
        payload["distance_stats"] = {k: float(v) for k, v in distance_stats.items()
                                     if not isinstance(v, (list, np.ndarray))}
    if calib is not None:
        payload["calibration"] = calib.to_dict()
    if extra:
        payload.update(extra)
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    return payload


def save_msd_npz(
    out_path: str,
    lags: np.ndarray,
    msd: np.ndarray,
    per_tile_msd: Optional[np.ndarray] = None,
    traj: Optional[np.ndarray] = None,
    per_tile_alpha: Optional[np.ndarray] = None,
    per_tile_K_alpha: Optional[np.ndarray] = None,
) -> None:
    """Write raw MSD arrays to NPZ for downstream re-analysis."""
    kwargs = {"lags": lags, "msd": msd}
    if per_tile_msd is not None:
        kwargs["per_tile_msd"] = per_tile_msd
    if traj is not None:
        kwargs["traj"] = traj
    if per_tile_alpha is not None:
        kwargs["per_tile_alpha"] = per_tile_alpha
    if per_tile_K_alpha is not None:
        kwargs["per_tile_K_alpha"] = per_tile_K_alpha
    np.savez_compressed(out_path, **kwargs)


# =============================================================================
# HIGH-LEVEL ONE-CALL WRAPPER
# =============================================================================

def run_msd_for_pair(
    conformations: Sequence[np.ndarray],
    pair: Tuple[int, int],
    label: str,
    out_dir: str,
    tile_size: int,
    pad: int,
    n_tiles: int,
    lag_min: int = 1,
    lag_max_frac: float = 1 / 3,
    fit_lag_min: int = 5,
    fit_lag_max_frac: float = 0.25,
    min_n_lags_for_fit: int = 6,
    calibration: Optional[Calibration] = None,
    save_npz: bool = True,
    file_prefix: str = "",
) -> Dict:
    """
    End-to-end: compute tiled two-point MSD for one pair, fit alpha, run
    saturation and distance-distribution diagnostics, save JSON (and NPZ).

    Returns a dict compatible with :func:`plot_msd_curves`, i.e. keys
    ``lags``, ``msd``, ``alpha``, ``D``, ``fit_range``. The same dict is
    persisted at ``<out_dir>/msd_<label>.json`` so downstream overlay
    plots can read back any subset of labels without recomputing.
    """
    os.makedirs(out_dir, exist_ok=True)
    idx_a, idx_b = pair
    t0 = time.time()
    result = compute_two_point_msd_tiled(
        conformations, idx_a=idx_a, idx_b=idx_b,
        tile_size=tile_size, pad=pad, n_tiles=n_tiles,
        lag_min=lag_min, lag_max_frac=lag_max_frac,
    )
    compute_secs = time.time() - t0

    lags = result["lags"]
    msd = result["msd"]
    fit = fit_msd_alpha(
        lags, msd,
        fit_lag_min=fit_lag_min,
        fit_lag_max_frac=fit_lag_max_frac,
        min_n_lags_for_fit=min_n_lags_for_fit,
    )
    sat_fit = fit_msd_saturation(lags, msd)

    # Per-tile fits: one (alpha, K_alpha) per tile.  Used for within-
    # replicate uncertainty bands and, at the pooling stage, for
    # between-condition statistics (restricted to one value per replicate).
    per_tile_fits = fit_msd_alpha_per_tile(
        lags, result.get("per_tile_msd"),
        fit_lag_min=fit_lag_min,
        fit_lag_max_frac=fit_lag_max_frac,
        min_n_lags_for_fit=min_n_lags_for_fit,
    )

    # Distance stats on the first tile's trajectory (enough for sanity checks).
    distance_stats: Optional[Dict[str, float]] = None
    try:
        chrom_size = conformations[0].shape[0]
        if chrom_size >= tile_size:
            a_abs0, b_abs0 = pad + idx_a, pad + idx_b
        else:
            a_abs0, b_abs0 = idx_a, idx_b
        traj0 = _pair_trajectory(conformations, a_abs0, b_abs0)
        distance_stats = pair_distance_distribution(traj0)
        # Drop histogram arrays from the summary dict; keep scalar stats.
        distance_stats = {k: v for k, v in distance_stats.items()
                          if k in ("mean_d", "median_d", "std_d")}
    except Exception as e:  # noqa: BLE001
        logger.debug(f"distance stats failed for {label}: {e}")

    # Strip ndarrays from per-tile fits before JSON serialisation; keep the
    # summary scalars there and the arrays in the companion NPZ.
    per_tile_scalars = {
        k: v for k, v in per_tile_fits.items()
        if k not in ("alpha_per_tile", "K_alpha_per_tile")
    }
    safe_label = label.replace("/", "_").replace(":", "_")
    fname_stem = f"{file_prefix}_msd_{safe_label}" if file_prefix else f"msd_{safe_label}"
    json_path = os.path.join(out_dir, f"{fname_stem}.json")
    save_msd_json(
        json_path, lags, msd, pair, label, fit,
        sat_fit=sat_fit, distance_stats=distance_stats, calib=calibration,
        extra={"n_tiles_used": int(result["n_tiles_used"]),
               "compute_seconds": compute_secs,
               "per_tile_fit": per_tile_scalars},
    )

    if save_npz:
        npz_path = os.path.join(out_dir, f"{fname_stem}.npz")
        save_msd_npz(
            npz_path, lags, msd,
            per_tile_msd=result.get("per_tile_msd"),
            per_tile_alpha=per_tile_fits.get("alpha_per_tile"),
            per_tile_K_alpha=per_tile_fits.get("K_alpha_per_tile"),
        )

    logger.info(
        f"  MSD({label}) pair=({idx_a},{idx_b}) sep={abs(idx_b-idx_a)} "
        f"monomers, n_tiles={result['n_tiles_used']}, "
        f"α={fit.get('alpha')}, K_α={fit.get('K_alpha')}, "
        f"compute={compute_secs:.2f}s"
    )

    return {
        "lags": lags,
        "msd": msd,
        "alpha": fit.get("alpha"),
        "K_alpha": fit.get("K_alpha"),
        "D": fit.get("D"),                   # legacy alias
        "alpha_stderr": fit.get("alpha_stderr"),
        "log_K_alpha_stderr": fit.get("log_K_alpha_stderr"),
        "fit_range": (fit.get("fit_lag_min"), fit.get("fit_lag_max"))
        if fit.get("alpha") is not None else None,
        "sat_fit": sat_fit,
        "pair": pair,
        "label": label,
        "json_path": json_path,
        "per_tile_fit": per_tile_fits,
    }
