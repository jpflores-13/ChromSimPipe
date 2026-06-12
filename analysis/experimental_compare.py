#!/usr/bin/env python
"""
Optional experimental Hi-C / Micro-C comparison helpers.

Loads a ``.mcool`` (or ``.cool``) file at the same resolution the simulation
uses, extracts the Sox2 locus, and provides convenience routines to compare
it with simulated data using the P(s) and APA modules.

If cooler is not installed or the mcool path doesn't exist, every function
here is a graceful no-op — callers should treat a None/empty return as
"no experimental comparison available".
"""

from __future__ import annotations

import json
import logging
import os
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# COOLER LOADING
# =============================================================================

def _cooler_available() -> bool:
    try:
        import cooler  # noqa: F401
        return True
    except ImportError:
        return False


def resolve_mcool_uri(path: str, resolution: int) -> str:
    """
    Return a proper ``::/resolutions/<bp>`` URI for an .mcool file.
    Passes .cool paths through unchanged.
    """
    if path.endswith(".mcool") and "::" not in path:
        return f"{path}::/resolutions/{int(resolution)}"
    return path


def load_region_matrix(
    mcool_path: str,
    chrom: str,
    start_bp: int,
    end_bp: int,
    resolution: int,
    balance: bool = True,
) -> Optional[np.ndarray]:
    """
    Fetch a dense (N, N) matrix for ``chrom:start-end`` at the given
    resolution. Returns ``None`` if cooler is missing or the file doesn't
    exist. NaNs (from balancing) are replaced with 0.

    ``N`` equals ``(end_bp - start_bp) // resolution`` — the caller can
    assume the simulated contact map and the experimental matrix share the
    same shape so long as ``resolution`` matches the simulation bin size.
    """
    if not _cooler_available():
        logger.warning("  cooler not installed — skipping experimental load")
        return None
    if not os.path.exists(mcool_path):
        logger.warning(f"  mcool not found: {mcool_path}")
        return None

    import cooler
    uri = resolve_mcool_uri(mcool_path, resolution)
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        logger.warning(f"  failed to open {uri}: {e}")
        return None

    region = f"{chrom}:{start_bp}-{end_bp}"
    try:
        mat = np.asarray(clr.matrix(balance=balance).fetch(region), dtype=float)
    except Exception as e:
        logger.warning(f"  failed to fetch {region} from {mcool_path}: {e}")
        return None

    mat[~np.isfinite(mat)] = 0.0

    # Pad/truncate to the expected size if cooler returns off-by-one
    N_expected = (end_bp - start_bp) // resolution
    N_actual = mat.shape[0]
    if N_actual != N_expected:
        logger.info(f"  mcool region returned {N_actual} bins, expected "
                    f"{N_expected} — trimming/padding")
        N = min(N_actual, N_expected)
        out = np.zeros((N_expected, N_expected), dtype=float)
        out[:N, :N] = mat[:N, :N]
        mat = out
    return mat


# =============================================================================
# SIDE-BY-SIDE PLOT: sim map vs exp map
# =============================================================================

def plot_sim_vs_exp_map(
    sim_map: np.ndarray,
    exp_map: np.ndarray,
    out_path: str,
    sim_label: str = "simulation",
    exp_label: str = "Hi-C",
    log_scale: bool = True,
    cmap: str = "fall",
) -> None:
    """
    Two-panel plot comparing simulated and experimental contact maps,
    both log-scaled to the same vmin/vmax.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        import cooltools  # noqa: F401
    except ImportError:
        if cmap == "fall":
            cmap = "Reds"

    N = min(sim_map.shape[0], exp_map.shape[0])
    sim = np.asarray(sim_map[:N, :N], dtype=float)
    exp = np.asarray(exp_map[:N, :N], dtype=float)
    if log_scale:
        eps_sim = max(np.nanmin(sim[sim > 0]) * 0.1 if np.any(sim > 0) else 1e-12, 1e-12)
        eps_exp = max(np.nanmin(exp[exp > 0]) * 0.1 if np.any(exp > 0) else 1e-12, 1e-12)
        sim = np.log10(sim + eps_sim)
        exp = np.log10(exp + eps_exp)

    vmin = float(np.nanpercentile(np.concatenate([sim.ravel(), exp.ravel()]), 2))
    vmax = float(np.nanpercentile(np.concatenate([sim.ravel(), exp.ravel()]), 99))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, mat, lbl in zip(axes, [sim, exp], [sim_label, exp_label]):
        im = ax.imshow(mat, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(lbl)
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02,
                 label="log10(contact)" if log_scale else "contact")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  sim-vs-exp contact map saved → {out_path}")


# =============================================================================
# SUMMARY METRICS
# =============================================================================

def compute_sim_exp_metrics(
    sim_map: np.ndarray,
    exp_map: np.ndarray,
    max_diag: Optional[int] = None,
) -> Dict[str, float]:
    """
    Pearson correlation per diagonal + overall, plus stratum-adjusted
    correlation and log-P(s) correlation. Mirrors
    ``contact_maps.compare_contact_maps`` but exposed as a standalone
    entry point so experimental_compare can be used on its own.
    """
    N = min(sim_map.shape[0], exp_map.shape[0])
    sim = sim_map[:N, :N].astype(float).copy()
    exp = exp_map[:N, :N].astype(float).copy()

    if max_diag is not None:
        r, c = np.ogrid[:N, :N]
        mask = np.abs(r - c) > max_diag
        sim[mask] = 0
        exp[mask] = 0

    sim_norm = sim / (np.nanmax(sim) + 1e-10)
    exp_norm = exp / (np.nanmax(exp) + 1e-10)

    diag_corrs: List[float] = []
    for d in range(1, min(N, max_diag or N)):
        sd = np.diagonal(sim_norm, offset=d)
        ed = np.diagonal(exp_norm, offset=d)
        if np.nanstd(sd) > 0 and np.nanstd(ed) > 0:
            diag_corrs.append(float(np.corrcoef(sd, ed)[0, 1]))

    upper = np.triu_indices(N, k=1)
    sim_u = sim_norm[upper]; exp_u = exp_norm[upper]
    ok = (sim_u > 0) | (exp_u > 0)
    overall = float(np.corrcoef(sim_u[ok], exp_u[ok])[0, 1]) if ok.sum() else 0.0

    weights = np.array([N - d for d in range(1, len(diag_corrs) + 1)], dtype=float)
    if weights.size:
        weights /= weights.sum()
        scc = float(np.sum(np.asarray(diag_corrs) * weights))
    else:
        scc = 0.0

    return {
        "overall_pearson": overall,
        "mean_diag_pearson": float(np.nanmean(diag_corrs)) if diag_corrs else 0.0,
        "stratum_adjusted_corr": scc,
        "n_diagonals_compared": int(len(diag_corrs)),
    }


def save_metrics_json(metrics: Dict[str, float], out_path: str) -> None:
    with open(out_path, "w") as f:
        json.dump(metrics, f, indent=2)
    logger.info(f"  sim-exp metrics saved → {out_path}")
