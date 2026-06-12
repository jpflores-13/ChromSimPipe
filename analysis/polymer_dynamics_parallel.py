#!/usr/bin/env python
"""
File-parallel, single-pass polymer-dynamics diagnostics.

A faithful numerical re-implementation of ``analysis/polymer_dynamics.py``
that collapses the seven sequential stream passes the orchestrator
currently performs (R_g + looped_fraction × n_pairs + dwell_times ×
n_pairs + pair_distance_dist × n_pairs) into **one** pass parallelised
across polychrom ``blocks_*.h5`` files — the same Hansen-lab worker
pattern already used by ``analysis.contact_maps._block_file_tile_contact_worker``.

Motivation
----------
On a pooled 3-replicate condition (~300 k frames × 28 tiles × 2000
beads) the existing sequential implementation spends most of Phase 2
in ``polymer_dynamics`` after MSD finishes: each of the 7 passes is a
single-threaded Python for-loop over all frames, streaming one frame
at a time through h5py. The 15 other cores allocated to the job sit
idle. For the ongoing job 9388554 this adds ~2 h per condition on
top of MSD, turning a would-fit 48 h sweep into a timeout.

This module addresses that by:

  1. Handing each ``blocks_*.h5`` file to a worker in a
     ``multiprocessing.Pool(n_jobs)``.
  2. Per block, each worker opens the file ONCE, walks all its frames
     once, and accumulates only the per-tile scalars needed by every
     downstream diagnostic: R_g, and pair separations for every
     requested pair.
  3. The parent concatenates per-block arrays in block-file order
     (preserving temporal ordering required by ``compute_dwell_times``)
     and runs the small, vectorised post-processing steps
     (histograms, contiguous-run detection, per-tile statistics).

Memory profile
--------------
Per worker:   ~1–3 MB scratch per frame + per-block output arrays that
              scale as n_tiles × n_frames_in_block (≤ a few MB per
              block). No full trajectory ever held.
Parent:       n_tiles × T × 8 B for R_g-per-tile and
              n_pairs × n_tiles × T × 8 B for pair distances.
              For T=300 k, n_tiles=28, n_pairs=2 this totals ~200 MB —
              three orders of magnitude below the current pooled
              contact-map footprint.

Numerical equivalence
---------------------
The per-tile R_g formula and the |r_b − r_a| distance formula are the
same as in ``polymer_dynamics.py`` (vectorised tile stacks, two-pass
R_g identity). Downstream post-processing functions
(``_contiguous_runs``, histogramming) are imported from there directly
so a regression in one file can't silently diverge the other.

Drop-in usage
-------------
``compute_dynamics_batch`` returns a single dict containing the same
fields that seven separate calls to the sequential API would produce::

    out["rg_timecourse"]                    # -> same as compute_rg_timecourse(...)
    out["per_pair"][label]["looped_fraction"]          # -> compute_looped_fraction(...)
    out["per_pair"][label]["dwell_times"]              # -> compute_dwell_times(...)
    out["per_pair"][label]["pair_distance_dist"]       # -> compute_pair_distance_distribution(...)

A thin ``run_polymer_dynamics_parallel`` wrapper mirrors the block that
``_run_msd_and_dynamics`` runs after MSD, so the orchestrator swap is
a one-line change.

Fallback
--------
If the input is not a ``ConformationStream`` with a populated
``block_files`` list (e.g. legacy ``conformations.h5`` format), the
module transparently falls through to the sequential implementation
in ``analysis.polymer_dynamics``, so nothing regresses on data it
can't parallelise.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from analysis.polymer_dynamics import (
    _contiguous_runs,
    compute_dwell_times as _seq_compute_dwell_times,
    compute_looped_fraction as _seq_compute_looped_fraction,
    compute_pair_distance_distribution as _seq_compute_pair_distance_distribution,
    compute_rg_timecourse as _seq_compute_rg_timecourse,
    plot_dwell_time_hist,
    plot_pair_distance_hist,
    plot_rg_timecourse,
)

logger = logging.getLogger(__name__)


# =============================================================================
# WORKER — one h5 block file at a time
# =============================================================================

def _per_block_dynamics_worker(args):
    """Hansen-lab file-parallel worker.

    Opens a single ``blocks_*.h5`` once, walks every frame inside it, and
    for each frame computes:

        (a) Per-tile R_g over a bead sub-range
        (b) Per-pair, per-tile |r_b - r_a| separation

    Returns
    -------
    rg        : (n_tiles, n_frames_in_block) float64
    pair_d    : (n_pairs, n_tiles, n_frames_in_block) float64
    n_frames  : int

    Returning numpy arrays (not Python lists) keeps the inter-process
    transfer cheap: ~1–3 MB per block even for a 5000-frame tiled
    trajectory with 28 tiles and 2 pairs.
    """
    import h5py
    import numpy as _np

    (block_file, tile_offsets, start_loc, end_loc, pair_indices) = args
    tile_offsets = _np.asarray(tile_offsets, dtype=int)
    pair_indices = _np.asarray(pair_indices, dtype=int)  # (n_pairs, 2)

    n_tiles = tile_offsets.size
    n_pairs = pair_indices.shape[0]
    N_sub = end_loc - start_loc

    with h5py.File(block_file, "r") as hf:
        keys = sorted(hf.keys(), key=lambda x: int(x) if x.isdigit() else 0)
        frame_keys = [k for k in keys
                      if isinstance(hf[k], h5py.Group) and "pos" in hf[k]]
        n_frames = len(frame_keys)

        if n_frames == 0:
            return (_np.empty((n_tiles, 0), dtype=_np.float64),
                    _np.empty((n_pairs, n_tiles, 0), dtype=_np.float64),
                    0)

        rg = _np.empty((n_tiles, n_frames), dtype=_np.float64)
        pair_d = _np.empty((n_pairs, n_tiles, n_frames), dtype=_np.float64)

        # Precompute the flat bead index set we need per frame to minimise
        # fancy-indexing overhead: the locus sub-range for every tile plus
        # the (idx_a, idx_b) ends for every requested pair, on every tile.
        # We read these slices directly per-frame rather than fancy-indexing
        # so the R_g einsum keeps its contiguous layout.
        for t_idx, key in enumerate(frame_keys):
            frame = _np.asarray(hf[key]["pos"], dtype=_np.float64)

            # ── R_g per tile (two-pass identity) ──────────────────────────
            # Stack n_tiles sub-slices of shape (N_sub, 3).
            subs = _np.stack(
                [frame[off + start_loc : off + end_loc, :] for off in tile_offsets],
                axis=0,
            )                                                    # (n_tiles, N_sub, 3)
            sum_x = subs.sum(axis=1)                             # (n_tiles, 3)
            sum_sq = _np.einsum("kij,kij->k", subs, subs)        # (n_tiles,)
            mean_sq = sum_sq / N_sub
            cm_sq = _np.einsum("kj,kj->k", sum_x, sum_x) / (N_sub * N_sub)
            rg[:, t_idx] = _np.sqrt(_np.maximum(mean_sq - cm_sq, 0.0))

            # ── Per-pair, per-tile separation ─────────────────────────────
            for p_idx in range(n_pairs):
                idx_a, idx_b = pair_indices[p_idx]
                # Stack the two beads across tiles: (n_tiles, 3) each.
                a_pos = frame[tile_offsets + idx_a]              # (n_tiles, 3)
                b_pos = frame[tile_offsets + idx_b]              # (n_tiles, 3)
                diff = b_pos - a_pos
                pair_d[p_idx, :, t_idx] = _np.linalg.norm(diff, axis=1)

    return rg, pair_d, n_frames


# =============================================================================
# PUBLIC API — batched, file-parallel entry point
# =============================================================================

def compute_dynamics_batch(
    conformations,
    *,
    pairs: Sequence[Tuple[int, int]],
    labels: Sequence[str],
    contact_radius: float,
    tile_size: int,
    pad: int,
    n_tiles: int,
    n_jobs: int,
    bead_range: Optional[Tuple[int, int]] = None,
    bins: int = 60,
) -> Dict[str, object]:
    """Single-pass, file-parallel R_g + per-pair dynamics.

    Parameters
    ----------
    conformations : ConformationStream (or compatible) with ``block_files``
        If ``block_files`` is not populated (legacy ``conformations.h5``,
        non-h5 stream, etc.) the function falls back to the sequential
        implementations in ``analysis.polymer_dynamics``.
    pairs, labels : parallel sequences
        Each pair is a (idx_a, idx_b) in *locus coordinates* (before
        adding ``tile_offset``). ``labels`` is the human label used to
        key the returned per-pair dict.
    contact_radius : float
        Used for looped-fraction and dwell-time binary trace.
    tile_size, pad, n_tiles : int
        Tiling parameters from ``configs.parameters.TILING``. Must be
        consistent with the trajectory (``chrom_size == tile_size * n_tiles``).
    n_jobs : int
        Number of workers. The dispatch unit is one ``blocks_*.h5``; a
        higher ``n_jobs`` than the number of block files just leaves some
        workers idle.
    bead_range : (start, end) in locus coordinates, optional
        Sub-range for R_g. Defaults to the whole locus.
    bins : int
        Number of histogram bins for the pair-distance distribution.

    Returns
    -------
    dict structured as::

        {
            "rg_timecourse": {"t": (T,) int, "rg": (T,) float,
                              "rg_per_tile": (n_tiles, T) float},
            "per_pair": {
                label: {
                    "looped_fraction": {...},
                    "dwell_times": {...},
                    "pair_distance_dist": {...},
                },
                ...
            },
            "T": int,
            "n_tiles_used": int,
            "source": "file-parallel" | "sequential-fallback",
        }

    Notes
    -----
    Temporal ordering is preserved by sorting block files
    lexicographically (matching ``load_conformations_h5``'s glob) and
    using ``Pool.map`` (FIFO) so dwell-time contiguous-run detection
    remains faithful to the on-disk frame order.
    """
    if len(pairs) != len(labels):
        raise ValueError(
            f"pairs (n={len(pairs)}) and labels (n={len(labels)}) must align.")

    # --- Check the trajectory shape is tiled as described -------------
    chrom_size = conformations[0].shape[0]
    if chrom_size != tile_size * n_tiles:
        raise ValueError(
            f"chrom_size={chrom_size} inconsistent with tile_size*n_tiles="
            f"{tile_size * n_tiles}; parallel path requires tiled data.")

    N_locus = tile_size - 2 * pad
    if bead_range is None:
        bead_range = (0, N_locus)
    start_loc, end_loc = bead_range
    if end_loc > N_locus:
        raise ValueError(
            f"bead_range end {end_loc} exceeds locus size {N_locus}.")

    block_files = getattr(conformations, "block_files", None)

    # --- Fallback: no block files, or single worker ------------------
    if not block_files or n_jobs <= 1:
        logger.info(
            "  [dynamics-parallel] falling back to sequential "
            "(block_files=%s, n_jobs=%d)",
            "present" if block_files else "none", n_jobs,
        )
        return _compute_dynamics_batch_sequential(
            conformations,
            pairs=pairs, labels=labels,
            contact_radius=contact_radius,
            tile_size=tile_size, pad=pad, n_tiles=n_tiles,
            bead_range=bead_range, bins=bins,
        )

    tile_offsets = (np.arange(n_tiles, dtype=int) * tile_size + pad).tolist()
    pair_list = [(int(a), int(b)) for a, b in pairs]

    args_list = [
        (bf, tile_offsets, int(start_loc), int(end_loc), pair_list)
        for bf in block_files
    ]

    n_workers = min(n_jobs, len(args_list))
    logger.info(
        "  [dynamics-parallel] %d block files × %d tiles × %d pair(s) "
        "dispatched to %d worker(s)",
        len(args_list), n_tiles, len(pair_list), n_workers,
    )

    # Use the "spawn" context so workers do not inherit the parent's
    # already-warm numpy / h5py state — keeps per-worker RSS bounded
    # over long Phase-2 sweeps. The file-parallel pattern makes the
    # spawn cost negligible (1 × per worker, not per task).
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers) as pool:
        # Pool.map preserves input order, which we need so that
        # concatenation along the frame axis matches on-disk order.
        results = pool.map(_per_block_dynamics_worker, args_list)

    # --- Concatenate in block-file order ------------------------------
    rg_parts = [r[0] for r in results if r[2] > 0]
    pd_parts = [r[1] for r in results if r[2] > 0]
    if not rg_parts:
        raise RuntimeError(
            "compute_dynamics_batch: no frames found in any block file "
            f"under {block_files[0]!r} ...")
    rg_per_tile = np.concatenate(rg_parts, axis=1)          # (n_tiles, T)
    pair_d = np.concatenate(pd_parts, axis=2)               # (n_pairs, n_tiles, T)
    T = rg_per_tile.shape[1]

    # --- R_g time-course output ---------------------------------------
    rg_out = {
        "t": np.arange(T, dtype=int),
        "rg": rg_per_tile.mean(axis=0),
        "rg_per_tile": rg_per_tile,
    }

    # --- Per-pair post-processing -------------------------------------
    per_pair: Dict[str, Dict[str, object]] = {}
    for p_idx, label in enumerate(labels):
        d_all = pair_d[p_idx]                               # (n_tiles, T)
        per_pair[label] = {
            "looped_fraction":
                _looped_fraction_from_distances(d_all, contact_radius),
            "dwell_times":
                _dwell_times_from_distances(d_all, contact_radius),
            "pair_distance_dist":
                _pair_distance_dist_from_distances(d_all, bins),
        }

    return {
        "rg_timecourse": rg_out,
        "per_pair": per_pair,
        "T": T,
        "n_tiles_used": int(n_tiles),
        "source": "file-parallel",
    }


# =============================================================================
# POST-PROCESSING HELPERS — vectorised, no stream access
# =============================================================================

def _looped_fraction_from_distances(
    d_all: np.ndarray, contact_radius: float,
) -> Dict[str, object]:
    """Same formula as ``compute_looped_fraction`` (tiled branch)."""
    fracs = (d_all < contact_radius).mean(axis=1)           # (n_tiles,)
    return {
        "looped_fraction": float(fracs.mean()),
        "looped_fraction_sem": (
            float(fracs.std(ddof=1) / np.sqrt(len(fracs)))
            if len(fracs) > 1 else 0.0
        ),
        "n_tiles_used": int(len(fracs)),
        "contact_radius": float(contact_radius),
        "per_tile_fraction": fracs.tolist(),
    }


def _dwell_times_from_distances(
    d_all: np.ndarray, contact_radius: float,
) -> Dict[str, object]:
    """Same contiguous-run collection as ``compute_dwell_times``."""
    n_tiles = d_all.shape[0]
    dwell_sets: List[np.ndarray] = [
        _contiguous_runs(d_all[k] < contact_radius) for k in range(n_tiles)
    ]
    dwells = (np.concatenate(dwell_sets) if dwell_sets
              else np.empty(0, dtype=int))
    if dwells.size == 0:
        return {
            "dwells_frames": [], "n_dwells": 0,
            "mean": None, "median": None, "max": None,
            "contact_radius": float(contact_radius),
        }
    return {
        "dwells_frames": dwells.tolist(),
        "n_dwells": int(dwells.size),
        "mean": float(dwells.mean()),
        "median": float(np.median(dwells)),
        "max": int(dwells.max()),
        "contact_radius": float(contact_radius),
    }


def _pair_distance_dist_from_distances(
    d_all: np.ndarray, bins: int,
) -> Dict[str, object]:
    """Same histogram as ``compute_pair_distance_distribution`` (tiled)."""
    d = d_all.ravel()
    counts, edges = np.histogram(d, bins=bins, density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return {
        "centers": centers,
        "density": counts,
        "mean_d": float(d.mean()),
        "median_d": float(np.median(d)),
        "std_d": float(d.std()),
    }


# =============================================================================
# SEQUENTIAL FALLBACK — reuses the existing, well-tested functions
# =============================================================================

def _compute_dynamics_batch_sequential(
    conformations,
    *,
    pairs: Sequence[Tuple[int, int]],
    labels: Sequence[str],
    contact_radius: float,
    tile_size: int,
    pad: int,
    n_tiles: int,
    bead_range: Tuple[int, int],
    bins: int,
) -> Dict[str, object]:
    """Last-resort path for non-polychrom trajectories.

    Calls the original sequential implementations one pair at a time.
    Slower (n_pairs × 3 + 1 stream passes vs 1 file-parallel pass) but
    bit-for-bit identical and requires no special stream metadata.
    """
    rg = _seq_compute_rg_timecourse(
        conformations, bead_range=bead_range,
        tile_size=tile_size, pad=pad, n_tiles=n_tiles,
    )
    per_pair: Dict[str, Dict[str, object]] = {}
    for (idx_a, idx_b), label in zip(pairs, labels):
        per_pair[label] = {
            "looped_fraction": _seq_compute_looped_fraction(
                conformations, idx_a, idx_b,
                contact_radius=contact_radius,
                tile_size=tile_size, pad=pad, n_tiles=n_tiles,
            ),
            "dwell_times": _seq_compute_dwell_times(
                conformations, idx_a, idx_b,
                contact_radius=contact_radius,
                tile_size=tile_size, pad=pad, n_tiles=n_tiles,
            ),
            "pair_distance_dist": _seq_compute_pair_distance_distribution(
                conformations, idx_a, idx_b,
                bins=bins,
                tile_size=tile_size, pad=pad, n_tiles=n_tiles,
            ),
        }
    return {
        "rg_timecourse": rg,
        "per_pair": per_pair,
        "T": int(rg["rg"].shape[0]),
        "n_tiles_used": int(n_tiles),
        "source": "sequential-fallback",
    }


# =============================================================================
# THIN CONVENIENCE WRAPPER — matches _run_msd_and_dynamics's plotting block
# =============================================================================

def run_polymer_dynamics_parallel(
    conformations,
    *,
    pairs: Sequence[Tuple[int, int]],
    labels: Sequence[str],
    out_dir: str,
    display_name: str,
    contact_radius: float,
    tile_size: int,
    pad: int,
    n_tiles: int,
    n_jobs: int,
    N_locus: Optional[int] = None,
    bins: int = 60,
    save_npz: bool = True,
) -> Dict[str, object]:
    """Drop-in for the polymer-dynamics block in ``_run_msd_and_dynamics``.

    Runs ``compute_dynamics_batch`` once, then emits the same on-disk
    artefacts the sequential code writes:

        rg_timecourse.png / rg_timecourse.npz
        loop_fractions.json
        dwell_times.png
        pair_distance.png

    Returns the full batched dict for callers that want to compose
    calibration / overlays downstream.
    """
    import json

    os.makedirs(out_dir, exist_ok=True)

    bead_range = (0, N_locus) if N_locus is not None else None
    out = compute_dynamics_batch(
        conformations,
        pairs=pairs, labels=labels,
        contact_radius=contact_radius,
        tile_size=tile_size, pad=pad, n_tiles=n_tiles,
        n_jobs=n_jobs,
        bead_range=bead_range,
        bins=bins,
    )

    rg = out["rg_timecourse"]

    # --- R_g time-course plot + npz -----------------------------------
    try:
        plot_rg_timecourse(
            {display_name: rg},
            os.path.join(out_dir, "rg_timecourse.png"),
            title=f"Radius of gyration: {display_name}",
        )
    except Exception as e:                                           # noqa: BLE001
        logger.warning("    [R_g] plot failed: %s", e)

    if save_npz:
        np.savez_compressed(
            os.path.join(out_dir, "rg_timecourse.npz"),
            t=rg["t"], rg=rg["rg"],
        )

    # --- Per-pair artefacts -------------------------------------------
    loop_fracs = []
    dwell_by_label = {}
    pair_dist_by_label = {}
    for label in labels:
        pp = out["per_pair"][label]
        loop_fracs.append({"label": label, **pp["looped_fraction"]})
        dwell_by_label[label] = pp["dwell_times"]
        pair_dist_by_label[label] = pp["pair_distance_dist"]

    if loop_fracs:
        with open(os.path.join(out_dir, "loop_fractions.json"), "w") as f:
            json.dump(loop_fracs, f, indent=2)

    if dwell_by_label:
        try:
            plot_dwell_time_hist(
                dwell_by_label,
                os.path.join(out_dir, "dwell_times.png"),
                title=f"Contact dwell times: {display_name}",
            )
        except Exception as e:                                       # noqa: BLE001
            logger.warning("    [dwell] plot failed: %s", e)

    if pair_dist_by_label:
        try:
            plot_pair_distance_hist(
                pair_dist_by_label,
                os.path.join(out_dir, "pair_distance.png"),
                title=f"Pair separation density: {display_name}",
            )
        except Exception as e:                                       # noqa: BLE001
            logger.warning("    [pair-dist] plot failed: %s", e)

    out["rg_mean"] = float(rg["rg"].mean())
    out["rg_std"] = float(rg["rg"].std())
    return out
