#!/usr/bin/env python
"""
Polymer-dynamics analyses that complement two-point MSD.

These are the "curated handful" of single-trajectory diagnostics most
commonly used to interpret cohesin-loop-extrusion simulations, in the
vocabulary of the Gabriele 2022, chromatin_dynamics, and Hansen-lab
AbsLoopQuant papers. All of them act on the same polychrom trajectory
(list of (N_beads, 3) conformation arrays) that the contact-map and
MSD modules already consume, so orchestration is trivial.

Analyses included
-----------------
1. Radius of gyration time-course
   R_g(t)^2 = (1/N) * sum_i |r_i(t) - r_cm(t)|^2
   Tracks how compact the locus is over simulation time. Strong loop
   extrusion should compact the polymer (lower R_g); loss of cohesin
   should expand it.

2. Looped-fraction for a pair (occupancy diagnostic)
   Fraction of frames where the pair separation is below a contact
   radius. Complements the APA / two-point MSD picture: "how often is
   the pair in contact?" vs. "what does its separation-vector dynamics
   look like?".

3. Contact dwell-time distribution for a pair
   Distribution of contiguous runs of frames where the pair is within
   contact radius. The mean dwell time is an in-silico analogue of the
   microscopy "looped residence time" (Gabriele 2022: 10-30 min median).

4. Single-point MSD for one monomer
   Classic single-particle MSD including center-of-mass drift. Compared
   with the two-point MSD of a pair, reveals how much of the apparent
   mobility comes from the two-body mode vs. bulk translation.

5. End-to-end distance distribution
   Histogram of |r(t)| across the trajectory, with mean, median, std.
   Useful as a one-slide descriptor of the pair's typical state.

Public API
----------
compute_rg_timecourse       : (T,) R_g per frame for a slice of beads.
plot_rg_timecourse          : overlay R_g across conditions.
compute_looped_fraction     : scalar fraction below contact radius.
compute_dwell_times         : list of contiguous-contact dwell lengths.
plot_dwell_time_hist        : overlay dwell-time histograms (log x).
compute_single_point_msd    : one-point MSD for one bead.
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

_CVD_PALETTE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
]


def _stream_selected_beads(
    conformations: Sequence[np.ndarray],
    bead_indices,
) -> np.ndarray:
    """Stream the trajectory once and return only the selected beads.

    For a ConformationStream over a tiled 70k-monomer simulation, replacing
    ``np.asarray(conformations)`` (which would allocate 168 GB for T=100k)
    with this helper keeps the full (T, N_chrom, 3) array off the parent
    heap — only the picked ``n_beads`` positions per frame are retained.

    Parameters
    ----------
    conformations : iterable of (N_chrom, 3) arrays with __len__
    bead_indices : array-like of int bead indices to retain per frame

    Returns
    -------
    ndarray of shape (T, n_beads, 3), dtype float64.
    Memory = T * n_beads * 24 B — e.g. 13 MB for T=100k, n_beads=56.
    """
    bead_indices = np.asarray(bead_indices, dtype=int)
    T = len(conformations)
    out = np.empty((T, bead_indices.size, 3), dtype=np.float64)
    for t, frame in enumerate(conformations):
        out[t] = frame[bead_indices]
    return out


# =============================================================================
# 1. RADIUS OF GYRATION TIME-COURSE
# =============================================================================

def compute_rg_timecourse(
    conformations: Sequence[np.ndarray],
    bead_range: Optional[Tuple[int, int]] = None,
    tile_size: Optional[int] = None,
    pad: Optional[int] = None,
    n_tiles: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    R_g(t) per frame, for the locus sub-range (handles tiling).

    Parameters
    ----------
    conformations : list of np.ndarray
        Each (N_beads, 3).
    bead_range : (start, end) in *locus coordinates*, optional
        Sub-range of beads to compute R_g over. Defaults to the whole locus.
    tile_size, pad, n_tiles : int, optional
        If all three are supplied AND the chromosome size equals
        ``tile_size * n_tiles``, R_g is averaged across tiles for the
        same locus sub-range.

    Returns
    -------
    dict with: ``t`` (frame index), ``rg`` (per-frame mean R_g),
    ``rg_per_tile`` (n_tiles, T) or None if not tiled.
    """
    if len(conformations) == 0:
        raise ValueError("Empty conformation list.")
    T = len(conformations)
    chrom_size = conformations[0].shape[0]

    is_tiled = (tile_size is not None and pad is not None and n_tiles is not None
                and chrom_size == tile_size * n_tiles)
    N_locus = tile_size - 2 * pad if is_tiled else chrom_size
    if bead_range is None:
        bead_range = (0, N_locus)
    start_loc, end_loc = bead_range
    if end_loc > N_locus:
        raise ValueError(f"bead_range end {end_loc} exceeds locus size {N_locus}")

    N_sub = end_loc - start_loc

    # Two-pass Rg identity:  Rg^2 = <|r|^2> - |<r>|^2
    # Streamed over frames so the (T, N_chrom, 3) trajectory is never
    # materialised (would be 168 GB for T=100k, N_chrom=70k).
    if is_tiled:
        tile_offsets = np.arange(n_tiles, dtype=int) * tile_size + pad
        per_tile = np.empty((n_tiles, T), dtype=np.float64)
        for t_idx, frame in enumerate(conformations):
            # Stack the n_tiles sub-slices for this frame: (n_tiles, N_sub, 3).
            # Total per-frame footprint ≈ n_tiles * N_sub * 24 B (1.3 MB for
            # 28 × 2000), released on next iteration.
            subs = np.stack(
                [frame[off + start_loc : off + end_loc, :] for off in tile_offsets],
                axis=0,
            )
            sum_x = subs.sum(axis=1)                        # (n_tiles, 3)
            sum_sq = np.einsum("kij,kij->k", subs, subs)    # (n_tiles,)
            mean_sq = sum_sq / N_sub
            cm_sq = np.einsum("kj,kj->k", sum_x, sum_x) / (N_sub * N_sub)
            per_tile[:, t_idx] = np.sqrt(np.maximum(mean_sq - cm_sq, 0.0))
        rg = per_tile.mean(axis=0)
        return {"t": np.arange(T, dtype=int), "rg": rg, "rg_per_tile": per_tile}

    # Non-tiled: stream per-frame scalars for the bead_range slice.
    rg = np.empty(T, dtype=np.float64)
    for t_idx, frame in enumerate(conformations):
        sub = frame[start_loc:end_loc, :]                   # (N_sub, 3)
        sum_x = sub.sum(axis=0)                             # (3,)
        sum_sq = float((sub * sub).sum())                   # scalar
        mean_sq = sum_sq / N_sub
        cm_sq = float(sum_x @ sum_x) / (N_sub * N_sub)
        rg[t_idx] = np.sqrt(max(mean_sq - cm_sq, 0.0))
    return {"t": np.arange(T, dtype=int), "rg": rg, "rg_per_tile": None}


def plot_rg_timecourse(
    curves: Dict[str, Dict[str, np.ndarray]],
    out_path: str,
    title: str = "Radius of gyration over simulation time",
    xlabel: str = "frame",
    ylabel: str = "R_g (monomer units)",
) -> None:
    """Overlay R_g(t) for multiple conditions on a linear axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for i, (label, d) in enumerate(curves.items()):
        color = _CVD_PALETTE[i % len(_CVD_PALETTE)]
        ax.plot(d["t"], d["rg"], "-", color=color, lw=1.4, alpha=0.85,
                label=f"{label} (mean {d['rg'].mean():.2f})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  R_g time-course saved: {out_path}")


# =============================================================================
# 2. LOOPED FRACTION (pair occupancy)
# =============================================================================

def _pair_distance_trajectory(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
) -> np.ndarray:
    """|r(t)| array for one pair (length T).

    Streams over frames, retaining only the two indexed positions per step,
    so the full (T, N_chrom, 3) trajectory is never materialised.
    """
    coords = _stream_selected_beads(conformations, [idx_a, idx_b])  # (T, 2, 3)
    diff = coords[:, 1, :] - coords[:, 0, :]                         # (T, 3)
    return np.linalg.norm(diff, axis=1).astype(np.float64)


def _pair_distance_trajectories_tiled(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    tile_size: int,
    pad: int,
    n_tiles: int,
) -> np.ndarray:
    """|r(t)| trajectories for the same (idx_a, idx_b) pair on every tile.

    Single streaming pass over the conformation iterable; output shape
    (n_tiles, T). Callers that previously looped ``n_tiles`` times calling
    ``_pair_distance_trajectory`` (which re-reads the HDF5 stream each
    call) should use this instead — one disk pass, ~27× faster.
    """
    tile_offsets = np.arange(n_tiles, dtype=int) * tile_size + pad
    bead_indices = np.empty(2 * n_tiles, dtype=int)
    bead_indices[0::2] = tile_offsets + idx_a
    bead_indices[1::2] = tile_offsets + idx_b
    coords = _stream_selected_beads(conformations, bead_indices)     # (T, 2*n_tiles, 3)
    T = coords.shape[0]
    coords = coords.reshape(T, n_tiles, 2, 3)
    diff = coords[:, :, 1, :] - coords[:, :, 0, :]                   # (T, n_tiles, 3)
    return np.linalg.norm(diff, axis=2).T.astype(np.float64)         # (n_tiles, T)


def compute_looped_fraction(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    contact_radius: float = 3.0,
    tile_size: Optional[int] = None,
    pad: Optional[int] = None,
    n_tiles: Optional[int] = None,
) -> Dict[str, float]:
    """
    Fraction of frames where the pair separation is within ``contact_radius``.

    If tiling is in use, every tile contributes one independent
    trajectory; the returned ``looped_fraction`` is the mean across tiles
    and ``looped_fraction_sem`` is the SEM across tiles.
    """
    if len(conformations) == 0:
        raise ValueError("Empty conformation list.")
    chrom_size = conformations[0].shape[0]
    is_tiled = (tile_size is not None and pad is not None and n_tiles is not None
                and chrom_size == tile_size * n_tiles)

    if not is_tiled:
        d = _pair_distance_trajectory(conformations, idx_a, idx_b)
        frac = float((d < contact_radius).mean())
        return {"looped_fraction": frac, "looped_fraction_sem": 0.0,
                "n_tiles_used": 1, "contact_radius": float(contact_radius)}

    # (n_tiles, T) in one stream pass — no per-tile disk re-read.
    d_all = _pair_distance_trajectories_tiled(
        conformations, idx_a, idx_b, tile_size, pad, n_tiles,
    )
    fracs = (d_all < contact_radius).mean(axis=1)                    # (n_tiles,)
    return {
        "looped_fraction": float(fracs.mean()),
        "looped_fraction_sem": float(fracs.std(ddof=1) / np.sqrt(len(fracs)))
        if len(fracs) > 1 else 0.0,
        "n_tiles_used": int(len(fracs)),
        "contact_radius": float(contact_radius),
        "per_tile_fraction": fracs.tolist(),
    }


# =============================================================================
# 3. CONTACT DWELL-TIME DISTRIBUTION
# =============================================================================

def compute_dwell_times(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    contact_radius: float = 3.0,
    tile_size: Optional[int] = None,
    pad: Optional[int] = None,
    n_tiles: Optional[int] = None,
) -> Dict[str, object]:
    """
    Distribution of contiguous run-lengths where the pair is in contact.

    The dwell times are collected from the binary contact trace (1 = within
    ``contact_radius``, 0 = outside). If tiling is in use, trajectories
    from all tiles are concatenated with a zero-break so that runs do not
    straddle tile boundaries.

    Returns a dict with dwell-time statistics and the raw dwell-time
    array.
    """
    chrom_size = conformations[0].shape[0]
    is_tiled = (tile_size is not None and pad is not None and n_tiles is not None
                and chrom_size == tile_size * n_tiles)

    dwell_sets: List[np.ndarray] = []
    if is_tiled:
        d_all = _pair_distance_trajectories_tiled(
            conformations, idx_a, idx_b, tile_size, pad, n_tiles,
        )                                                            # (n_tiles, T)
        for k in range(n_tiles):
            dwell_sets.append(_contiguous_runs(d_all[k] < contact_radius))
    else:
        d = _pair_distance_trajectory(conformations, idx_a, idx_b)
        dwell_sets.append(_contiguous_runs(d < contact_radius))

    dwells = np.concatenate(dwell_sets) if dwell_sets else np.empty(0, dtype=int)
    if dwells.size == 0:
        return {"dwells_frames": [], "n_dwells": 0, "mean": None,
                "median": None, "max": None, "contact_radius": float(contact_radius)}
    return {
        "dwells_frames": dwells.tolist(),
        "n_dwells": int(dwells.size),
        "mean": float(dwells.mean()),
        "median": float(np.median(dwells)),
        "max": int(dwells.max()),
        "contact_radius": float(contact_radius),
    }


def _contiguous_runs(mask: np.ndarray) -> np.ndarray:
    """Return the lengths of contiguous True runs in a 1D boolean array."""
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0:
        return np.empty(0, dtype=int)
    # Pad with False so the diff trick catches leading/trailing runs.
    padded = np.concatenate([[False], mask, [False]])
    edges = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return ends - starts


def plot_dwell_time_hist(
    per_label: Dict[str, Dict[str, object]],
    out_path: str,
    title: str = "Contact dwell-time distribution",
    bins: int = 40,
    log_x: bool = True,
) -> None:
    """Overlay dwell-time histograms (CVD-safe palette)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    all_dwells = []
    for d in per_label.values():
        all_dwells.extend(d.get("dwells_frames", []))
    if not all_dwells:
        logger.warning("  No dwell times to plot.")
        return
    arr = np.asarray(all_dwells, dtype=float)
    if log_x:
        edges = np.logspace(np.log10(max(arr.min(), 1)),
                             np.log10(arr.max()), bins + 1)
    else:
        edges = np.linspace(arr.min(), arr.max(), bins + 1)

    for i, (label, d) in enumerate(per_label.items()):
        dwells = np.asarray(d.get("dwells_frames", []), dtype=float)
        if dwells.size == 0:
            continue
        color = _CVD_PALETTE[i % len(_CVD_PALETTE)]
        hist, _ = np.histogram(dwells, bins=edges, density=True)
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax.plot(centers, hist, "-", lw=1.6, color=color,
                label=f"{label} (n={dwells.size}, mean {dwells.mean():.1f})")
    if log_x:
        ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dwell time (frames)")
    ax.set_ylabel("density")
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  Dwell-time histogram saved: {out_path}")


# =============================================================================
# 4. SINGLE-POINT MSD
# =============================================================================

def compute_single_point_msd(
    conformations: Sequence[np.ndarray],
    idx: int,
    lag_min: int = 1,
    lag_max: Optional[int] = None,
    lag_max_frac: float = 1 / 3,
    tile_size: Optional[int] = None,
    pad: Optional[int] = None,
    n_tiles: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    One-point MSD of a single monomer::

        MSD_1pt(tau) = < |r(t+tau) - r(t)|^2 >_t

    Unlike two-point MSD, this includes center-of-mass drift. Useful as a
    consistency check and for direct comparison with single-molecule
    tracking (where each tag is reported independently).
    """
    if len(conformations) == 0:
        raise ValueError("Empty conformation list.")
    T = len(conformations)
    chrom_size = conformations[0].shape[0]
    is_tiled = (tile_size is not None and pad is not None and n_tiles is not None
                and chrom_size == tile_size * n_tiles)

    def _msd_one(traj: np.ndarray, lmax: int, lmin: int = 1):
        lags = np.arange(lmin, lmax + 1, dtype=int)
        msd = np.empty(lags.shape, dtype=np.float64)
        for k, tau in enumerate(lags):
            diff = traj[tau:] - traj[:-tau]
            msd[k] = np.einsum("ij,ij->i", diff, diff).mean()
        return lags, msd

    if lag_max is None:
        lag_max = max(1, int(T * lag_max_frac))

    if not is_tiled:
        coords = _stream_selected_beads(conformations, [idx])            # (T, 1, 3)
        traj = coords[:, 0, :]
        lags, msd = _msd_one(traj, lag_max, lag_min)
        return {"lags": lags, "msd": msd, "n_tiles_used": 1}

    # One streaming pass: stack bead ``idx`` from every tile simultaneously.
    tile_offsets = np.arange(n_tiles, dtype=int) * tile_size + pad
    bead_indices = tile_offsets + idx                                    # (n_tiles,)
    coords = _stream_selected_beads(conformations, bead_indices)         # (T, n_tiles, 3)

    per_tile = []
    lags_ref: Optional[np.ndarray] = None
    for tile in range(n_tiles):
        traj = coords[:, tile, :]                                        # (T, 3)
        lags, msd = _msd_one(traj, lag_max, lag_min)
        per_tile.append(msd)
        lags_ref = lags if lags_ref is None else lags_ref
    arr = np.vstack(per_tile)
    return {"lags": lags_ref, "msd": arr.mean(axis=0),
            "n_tiles_used": int(n_tiles), "per_tile_msd": arr}


# =============================================================================
# 5. END-TO-END / PAIR DISTANCE DISTRIBUTION
# =============================================================================

def compute_pair_distance_distribution(
    conformations: Sequence[np.ndarray],
    idx_a: int,
    idx_b: int,
    bins: int = 60,
    tile_size: Optional[int] = None,
    pad: Optional[int] = None,
    n_tiles: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    Steady-state histogram of the pair separation |r(t)|, pooled across
    tiles if the chromosome is tiled.
    """
    chrom_size = conformations[0].shape[0]
    is_tiled = (tile_size is not None and pad is not None and n_tiles is not None
                and chrom_size == tile_size * n_tiles)

    if is_tiled:
        # (n_tiles, T) in one stream pass; flatten to pool across tiles.
        d_all = _pair_distance_trajectories_tiled(
            conformations, idx_a, idx_b, tile_size, pad, n_tiles,
        )
        d = d_all.ravel()
    else:
        d = _pair_distance_trajectory(conformations, idx_a, idx_b)
    counts, edges = np.histogram(d, bins=bins, density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return {
        "centers": centers,
        "density": counts,
        "mean_d": float(d.mean()),
        "median_d": float(np.median(d)),
        "std_d": float(d.std()),
    }


def plot_pair_distance_hist(
    per_label: Dict[str, Dict[str, np.ndarray]],
    out_path: str,
    title: str = "Pair separation distribution (steady state)",
    xlabel: str = "|r| (monomer units)",
) -> None:
    """Overlay pair-distance histograms across conditions/pairs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for i, (label, d) in enumerate(per_label.items()):
        color = _CVD_PALETTE[i % len(_CVD_PALETTE)]
        ax.plot(d["centers"], d["density"], "-", lw=1.6, color=color,
                label=f"{label} (mean {d['mean_d']:.2f})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  Pair-distance histogram saved: {out_path}")
