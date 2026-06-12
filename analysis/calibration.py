#!/usr/bin/env python
"""
Physical-units calibration for the cohesin-simulation project.

Our simulations run in dimensionless "monomer" units for distance and
"saved frame" units for time. To talk to microscopy experiments we need
real numbers: how many nanometres is one monomer, how many seconds is
one frame. This module produces that mapping.

Two anchor strategies, selectable at runtime
--------------------------------------------
* ``"hic"`` — match the simulation P(s) curve to an experimental Hi-C
  P(s) at a chosen reference separation. Works whenever you already have
  a Bonev mcool sitting on disk; no microscopy required.

* ``"msd"`` — Gabriele/Fbn2 formulas: match the steady-state plateau and
  short-time diffusivity of the two-point MSD between simulation and
  experiment to solve for dx_sim (nm) and dt_sim (s). Needs a microscopy
  MSD curve for the same genomic separation.

Either returns a :class:`Calibration` dataclass with two scalars and a
``source`` tag so plots can be labelled truthfully.

Users pick the anchor at runtime (``--calibrate-with hic|msd``); the
default is ``hic`` because it uses data we always have.

Notes on the formulas
---------------------
Hi-C match: under the assumption that simulation and experiment share
the same P(s) shape but differ in amplitude, the ratio of P(s) at a
reference distance s_ref gives the renormalisation of the contact
probability. For a fully-calibrated monomer size one also needs a
persistence-length estimate (~50 nm for 1 kb chromatin in mESCs, per
Gabriele 2022).

MSD match (Fbn2 style):
    dx_sim = sqrt(J_expt / J_sim) * dx_expt
    dt_sim = ( (G_sim * dx_sim^2) / (G_expt * dx_expt^2) )^2 * dt_expt
where J is the steady-state plateau and G is the short-time diffusivity.

Public API
----------
Calibration            : dataclass with nm_per_monomer, sec_per_frame, source, metadata.
get_calibration        : dispatcher (anchor = "hic" | "msd").
calibrate_from_hic     : match sim P(s) to experimental P(s) at s_ref.
calibrate_from_msd     : Fbn2-style MSD match.
apply_calibration_msd  : convert (lag_frames, msd_monomer2) to (sec, nm^2).
load_experimental_msd_csv : parse a simple experimental MSD CSV.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field, asdict
from typing import Dict, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class Calibration:
    """Physical-units mapping between simulation and real-world scales."""
    nm_per_monomer: float
    sec_per_frame: float
    source: str                         # "hic" | "msd" | "assumed"
    metadata: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


# =============================================================================
# Hi-C P(s) anchor
# =============================================================================

def calibrate_from_hic(
    sim_distances_bp: np.ndarray,
    sim_ps: np.ndarray,
    expt_distances_bp: np.ndarray,
    expt_ps: np.ndarray,
    s_ref_bp: int = 100_000,
    bp_per_monomer: int = 1000,
    persistence_nm: float = 50.0,
    sec_per_frame: float = 1.0,
) -> Calibration:
    """
    Anchor nm-per-monomer to the experimental P(s) amplitude at ``s_ref_bp``.

    Logic
    -----
    1. Interpolate log P(s) at ``s_ref_bp`` for sim and for experiment.
    2. The amplitude ratio ``A = P_expt(s_ref) / P_sim(s_ref)`` is a scalar
       that absorbs the sim contact-radius convention (sim P(s) counts
       monomer-pairs within contact_radius; experiment counts reads).
       It does *not* change the nm/bp conversion.
    3. nm-per-monomer is set from the polymer-physics estimate:
           nm_per_mon = persistence_nm * sqrt(bp_per_monomer / 1000)
       which reduces to ``persistence_nm`` at 1 kb/monomer. The amplitude
       ratio is recorded in metadata so downstream plots can show
       un-normalised contact frequencies.
    4. ``sec_per_frame`` is taken as-is: this anchor does not constrain
       time. Callers who also want a time calibration should pair this
       with :func:`calibrate_from_msd`.

    Parameters
    ----------
    sim_distances_bp, sim_ps : np.ndarray
        Simulated P(s) curve (from ``analysis.ps_curve.compute_ps_curve``).
    expt_distances_bp, expt_ps : np.ndarray
        Experimental P(s) curve (from a cooler file via
        ``compute_ps_from_cooler``).
    s_ref_bp : int
        Genomic separation where the amplitude ratio is taken. 100 kb is a
        robust default: well past the short-range bumps, well before the
        long-range noise floor.
    bp_per_monomer : int
        Base pairs per monomer in the simulation.
    persistence_nm : float
        Estimated persistence length of the coarse-grained polymer (nm).

    Returns
    -------
    Calibration
    """
    sim_ratio = _interp_log_y_at_x(sim_distances_bp, sim_ps, s_ref_bp)
    exp_ratio = _interp_log_y_at_x(expt_distances_bp, expt_ps, s_ref_bp)

    amplitude_ratio = np.nan
    if (sim_ratio is not None and exp_ratio is not None
            and sim_ratio > 0 and np.isfinite(sim_ratio)):
        amplitude_ratio = float(exp_ratio / sim_ratio)

    nm_per_mon = persistence_nm * np.sqrt(max(bp_per_monomer, 1) / 1000.0)
    return Calibration(
        nm_per_monomer=float(nm_per_mon),
        sec_per_frame=float(sec_per_frame),
        source="hic",
        metadata={
            "s_ref_bp": float(s_ref_bp),
            "bp_per_monomer": float(bp_per_monomer),
            "persistence_nm": float(persistence_nm),
            "amplitude_ratio_expt_over_sim": float(amplitude_ratio),
        },
    )


def _interp_log_y_at_x(x: np.ndarray, y: np.ndarray, x_target: float) -> Optional[float]:
    """Log-log interpolation: return y(x_target) from (x, y) arrays."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = (x > 0) & (y > 0) & np.isfinite(y)
    if valid.sum() < 2:
        return None
    lx = np.log10(x[valid])
    ly = np.log10(y[valid])
    if x_target <= 0:
        return None
    lxt = np.log10(x_target)
    if lxt < lx.min() or lxt > lx.max():
        return None
    return float(10.0 ** np.interp(lxt, lx, ly))


# =============================================================================
# MSD anchor (Gabriele / Fbn2 formulas)
# =============================================================================

def calibrate_from_msd(
    sim_lags_frames: np.ndarray,
    sim_msd_monomer2: np.ndarray,
    expt_lags_sec: np.ndarray,
    expt_msd_um2: np.ndarray,
    fit_lag_min_sec: float = 1.0,
    fit_lag_max_sec: float = 30.0,
) -> Calibration:
    """
    Match simulation MSD to experimental MSD using the Fbn2 calibration.

    For short-time diffusivity G and steady-state plateau J, define:

        dx_sim  = sqrt(J_expt / J_sim) * dx_expt            (nm/monomer, with dx_expt = 1 um)
        dt_sim  = ((G_sim*dx_sim^2) / (G_expt*dx_expt^2))^2 * dt_expt

    G is estimated by a log-log linear fit in the short-time window
    ``[fit_lag_min_sec, fit_lag_max_sec]``; J is taken as ``max(MSD)/2``
    (a robust conservative estimate; Bayesian fits from
    ``bayesmsd.TwoLocusRouseFit`` are a more refined alternative when the
    dependency is available).

    All lengths are returned in nm.
    """
    G_expt = _short_time_prefactor(expt_lags_sec, expt_msd_um2,
                                    fit_lag_min_sec, fit_lag_max_sec)
    # Use same window in frame units for sim (caller can rescale later)
    G_sim = _short_time_prefactor(sim_lags_frames, sim_msd_monomer2,
                                    fit_lag_min_sec, fit_lag_max_sec)

    sim_valid = (sim_lags_frames > 0) & (sim_msd_monomer2 > 0) & np.isfinite(sim_msd_monomer2)
    expt_valid = (expt_lags_sec > 0) & (expt_msd_um2 > 0) & np.isfinite(expt_msd_um2)

    if sim_valid.sum() < 3 or expt_valid.sum() < 3:
        raise ValueError("Not enough valid MSD samples to calibrate.")

    J_sim = float(np.nanmax(sim_msd_monomer2[sim_valid])) / 2.0
    J_expt = float(np.nanmax(expt_msd_um2[expt_valid])) / 2.0

    if J_sim <= 0 or J_expt <= 0:
        raise ValueError("Degenerate plateau: J_sim or J_expt is zero.")

    dx_sim_um = float(np.sqrt(J_expt / J_sim))           # um per monomer
    if G_sim is None or G_expt is None or G_sim <= 0 or G_expt <= 0:
        raise ValueError("Could not fit short-time diffusivity; "
                         "check your lag-window choice.")
    dt_sim_sec = float(((G_sim * dx_sim_um ** 2) /
                        (G_expt * 1.0 ** 2)) ** 2)
    return Calibration(
        nm_per_monomer=float(dx_sim_um * 1000.0),
        sec_per_frame=float(dt_sim_sec),
        source="msd",
        metadata={
            "J_sim_monomer2": J_sim,
            "J_expt_um2": J_expt,
            "G_sim_monomer2_per_frame": float(G_sim),
            "G_expt_um2_per_sec": float(G_expt),
            "fit_lag_min_sec": float(fit_lag_min_sec),
            "fit_lag_max_sec": float(fit_lag_max_sec),
        },
    )


def _short_time_prefactor(x: np.ndarray, y: np.ndarray,
                            x_min: float, x_max: float) -> Optional[float]:
    """Log-log fit: return 10**intercept (the prefactor) over [x_min, x_max]."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = (x >= x_min) & (x <= x_max) & (y > 0) & np.isfinite(y)
    if mask.sum() < 3:
        return None
    slope, intercept = np.polyfit(np.log10(x[mask]), np.log10(y[mask]), deg=1)
    return float(10.0 ** intercept)


# =============================================================================
# Runtime dispatcher
# =============================================================================

def get_calibration(anchor: str = "hic", **kwargs) -> Calibration:
    """
    Pick a calibration strategy based on the runtime anchor choice.

    Parameters
    ----------
    anchor : {"hic", "msd", "assumed"}
        ``"hic"``: forwards to :func:`calibrate_from_hic`. Requires
        ``sim_distances_bp, sim_ps, expt_distances_bp, expt_ps, s_ref_bp``.

        ``"msd"``: forwards to :func:`calibrate_from_msd`. Requires
        ``sim_lags_frames, sim_msd_monomer2, expt_lags_sec, expt_msd_um2``.

        ``"assumed"``: skip fitting. Returns a bare :class:`Calibration`
        with ``nm_per_monomer`` and ``sec_per_frame`` from ``kwargs``.

    Any other keywords are forwarded to the matching helper.
    """
    anchor = (anchor or "hic").lower()
    if anchor == "hic":
        return calibrate_from_hic(**kwargs)
    if anchor == "msd":
        return calibrate_from_msd(**kwargs)
    if anchor == "assumed":
        return Calibration(
            nm_per_monomer=float(kwargs.get("nm_per_monomer", 50.0)),
            sec_per_frame=float(kwargs.get("sec_per_frame", 1.0)),
            source="assumed",
            metadata={},
        )
    raise ValueError(f"Unknown calibration anchor {anchor!r}; "
                     f"expected 'hic', 'msd', or 'assumed'.")


# =============================================================================
# Unit conversions
# =============================================================================

def apply_calibration_msd(
    lags_frames: np.ndarray,
    msd_monomer2: np.ndarray,
    calib: Calibration,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert an MSD curve from (frame, monomer^2) to (second, nm^2)."""
    lags_sec = np.asarray(lags_frames, dtype=float) * calib.sec_per_frame
    msd_nm2 = np.asarray(msd_monomer2, dtype=float) * (calib.nm_per_monomer ** 2)
    return lags_sec, msd_nm2


# =============================================================================
# I/O helpers
# =============================================================================

def load_experimental_msd_csv(
    csv_path: str,
    col_time: str = "dt_s",
    col_msd: str = "msd_um2",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read a two-column CSV (time_in_seconds, MSD_in_um^2) into numpy arrays.

    Accepts either quoted ``"dt_s","msd_um2"`` headers or the exact column
    names passed via ``col_time`` / ``col_msd``. Tabs and commas both fine.
    """
    with open(csv_path) as f:
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        reader = csv.DictReader(f, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError(f"{csv_path}: no header row found")
        # Case-insensitive column matching
        mapping = {k.strip().lower(): k for k in reader.fieldnames}
        if col_time.lower() not in mapping or col_msd.lower() not in mapping:
            raise ValueError(
                f"{csv_path}: expected columns {col_time!r} and {col_msd!r}, "
                f"found {reader.fieldnames!r}"
            )
        ts, ms = [], []
        for row in reader:
            try:
                ts.append(float(row[mapping[col_time.lower()]]))
                ms.append(float(row[mapping[col_msd.lower()]]))
            except (ValueError, KeyError):
                continue
    return np.asarray(ts, dtype=float), np.asarray(ms, dtype=float)
