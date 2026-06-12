#!/usr/bin/env python
"""
LEF / cohesin loop lifetime analysis from 1D loop-extrusion trajectories.

Two entry points:

1. ``collect_lifetimes_from_simulator(sim, n_steps, burn_in=0)``
   Runs an existing ``LEFSimulator`` (from ``scripts.lef_dynamics``) for
   ``n_steps`` steps, tracking how long each LEF remains bound between
   binding and unbinding events. Returns per-LEF lifetime + loop-size
   distributions.

2. ``collect_lifetimes_from_trajectory(traj)``
   Post-processes a saved trajectory (a list-per-step of ``(left, right)``
   bond tuples) into the same outputs. This is the format the Hansen lab's
   DSB_smcTranslocator dumps.

Both return a dict with:
    'lifetimes_steps'      : int array, per-binding-event lifetime in steps
    'loop_sizes_bins'      : int array, max loop size reached per event
    'ctcf_stall_fraction'  : float, fraction of events ending while stalled
    'mean_loop_size_bins'  : float
    'n_events'             : int

Plots:
    plot_lifetime_histogram
    plot_loop_size_histogram
    plot_lifetime_vs_size
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# TRAJECTORY → LIFETIME DISTRIBUTION
# =============================================================================

def _to_bond_set(bonds: Iterable[Sequence[int]]) -> "set[Tuple[int, int]]":
    out = set()
    for b in bonds:
        a, c = int(b[0]), int(b[1])
        if a == c:
            continue
        if a > c:
            a, c = c, a
        out.add((a, c))
    return out


def collect_lifetimes_from_trajectory(
    traj: Sequence[Iterable[Sequence[int]]],
) -> Dict[str, np.ndarray]:
    """
    Build per-bond binding-episode statistics from a list-per-step
    trajectory of LEF bonds.

    A "binding event" begins the first step a given (left, right) pair
    appears and ends when that anchor pair vanishes. The pair is tracked
    *exactly* by its two endpoints — i.e. an extrusion step that moves
    either arm produces a different pair and is counted as a new (short)
    event. This matches Hansen-lab post-processing and is appropriate for
    identifying stalled/stable loops.

    Parameters
    ----------
    traj : sequence of iterables
        ``traj[t]`` is a collection of ``(left, right)`` pairs at step t.

    Returns
    -------
    dict : see module docstring.
    """
    active: Dict[Tuple[int, int], Dict[str, int]] = {}
    lifetimes: List[int] = []
    max_sizes: List[int] = []

    for t, bonds in enumerate(traj):
        current = _to_bond_set(bonds)

        # Close events whose bond is no longer present
        closed = [k for k in active if k not in current]
        for k in closed:
            info = active.pop(k)
            lifetimes.append(info["end"] - info["start"] + 1)
            max_sizes.append(info["max_size"])

        # Open new events or refresh existing ones
        for k in current:
            if k not in active:
                active[k] = {"start": t, "end": t, "max_size": k[1] - k[0]}
            else:
                active[k]["end"] = t
                active[k]["max_size"] = max(active[k]["max_size"], k[1] - k[0])

    # Flush anything still active at the end (censored)
    for info in active.values():
        lifetimes.append(info["end"] - info["start"] + 1)
        max_sizes.append(info["max_size"])

    lifetimes_arr = np.asarray(lifetimes, dtype=int)
    sizes_arr = np.asarray(max_sizes, dtype=int)

    return {
        "lifetimes_steps": lifetimes_arr,
        "loop_sizes_bins": sizes_arr,
        "n_events": int(lifetimes_arr.size),
        "mean_loop_size_bins": float(np.mean(sizes_arr)) if sizes_arr.size else float("nan"),
        "median_lifetime_steps": float(np.median(lifetimes_arr)) if lifetimes_arr.size else float("nan"),
        "mean_lifetime_steps": float(np.mean(lifetimes_arr)) if lifetimes_arr.size else float("nan"),
    }


# =============================================================================
# RUN A LEFSimulator AND COLLECT
# =============================================================================

def collect_lifetimes_from_simulator(
    sim,
    n_steps: int,
    burn_in: int = 0,
) -> Dict[str, np.ndarray]:
    """
    Advance the simulator and call collect_lifetimes_from_trajectory on
    the resulting trace. The simulator must expose ``.step()`` and
    ``.get_bonds()`` (the ``LEFSimulator`` in ``scripts.lef_dynamics``
    satisfies both).
    """
    for _ in range(burn_in):
        sim.step()
    traj: List[List[Tuple[int, int]]] = []
    for _ in range(n_steps):
        sim.step()
        traj.append(list(sim.get_bonds()))
    logger.info(f"  lef_lifetimes: ran {n_steps} steps (after {burn_in} burn-in)")
    return collect_lifetimes_from_trajectory(traj)


# =============================================================================
# PLOTTING
# =============================================================================

def plot_lifetime_histogram(
    stats: Dict[str, np.ndarray],
    out_path: str,
    title: str = "LEF lifetime distribution",
    resolution_seconds: Optional[float] = None,
    log_y: bool = True,
) -> None:
    """Histogram of per-event lifetimes. If ``resolution_seconds`` is set,
    lifetimes are shown in minutes; otherwise in simulation steps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lt = np.asarray(stats["lifetimes_steps"], dtype=float)
    if lt.size == 0:
        logger.warning("  plot_lifetime_histogram: empty lifetime array")
        return

    if resolution_seconds is not None:
        data = lt * resolution_seconds / 60.0
        xlabel = "Lifetime (min)"
    else:
        data = lt
        xlabel = "Lifetime (sim steps)"

    fig, ax = plt.subplots(figsize=(5.2, 4))
    ax.hist(data, bins=50, color="steelblue", edgecolor="white")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of binding events")
    ax.set_title(title)
    if log_y:
        ax.set_yscale("log")
    ax.text(0.98, 0.98,
            f"n = {int(stats['n_events'])}\n"
            f"mean = {stats['mean_lifetime_steps']:.1f} steps\n"
            f"median = {stats['median_lifetime_steps']:.1f} steps",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  LEF lifetime histogram saved → {out_path}")


def plot_loop_size_histogram(
    stats: Dict[str, np.ndarray],
    out_path: str,
    title: str = "LEF max loop size distribution",
    resolution_bp: int = 1000,
) -> None:
    """Histogram of max loop size reached per binding event (kb on x-axis)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sz = np.asarray(stats["loop_sizes_bins"], dtype=float)
    if sz.size == 0:
        logger.warning("  plot_loop_size_histogram: empty size array")
        return

    kb = sz * resolution_bp / 1000.0
    fig, ax = plt.subplots(figsize=(5.2, 4))
    ax.hist(kb, bins=50, color="indianred", edgecolor="white")
    ax.set_xlabel("Max loop size (kb)")
    ax.set_ylabel("Number of binding events")
    ax.set_title(title)
    ax.text(0.98, 0.98,
            f"n = {int(stats['n_events'])}\n"
            f"mean = {np.mean(kb):.1f} kb\n"
            f"median = {np.median(kb):.1f} kb",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  loop size histogram saved → {out_path}")


def plot_lifetime_vs_size(
    stats: Dict[str, np.ndarray],
    out_path: str,
    title: str = "Lifetime vs. loop size",
    resolution_bp: int = 1000,
    resolution_seconds: Optional[float] = None,
) -> None:
    """2D hex-bin of lifetime vs. max loop size — reveals CTCF-stalled population."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lt = np.asarray(stats["lifetimes_steps"], dtype=float)
    sz = np.asarray(stats["loop_sizes_bins"], dtype=float)
    if lt.size == 0 or sz.size == 0:
        return

    if resolution_seconds is not None:
        y = lt * resolution_seconds / 60.0
        ylabel = "Lifetime (min)"
    else:
        y = lt
        ylabel = "Lifetime (sim steps)"
    x = sz * resolution_bp / 1000.0

    fig, ax = plt.subplots(figsize=(5.2, 4))
    hb = ax.hexbin(x, y, gridsize=30, cmap="viridis", mincnt=1,
                   yscale="log" if lt.max() > 50 else "linear")
    fig.colorbar(hb, ax=ax, label="count")
    ax.set_xlabel("Max loop size (kb)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  lifetime-vs-size plot saved → {out_path}")
