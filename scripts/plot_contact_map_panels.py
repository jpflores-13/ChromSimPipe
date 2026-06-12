#!/usr/bin/env python
"""Generate contact-map figure panels from results/analysis/.

Reads the prefixed ``.npy`` contact maps written by run_analysis_all.py and
emits three families of figures into ``results/figures/`` (or wherever you
point ``--output-dir``):

  1. **Per-condition per-rep grid** — one subplot per replicate of one
     condition, drawn as log-scale heatmaps on a shared colour scale.
     File: ``{condition}_contact_maps_per_rep.png``

  2. **Per-condition pooled (merged-reps)** — single contact-map heatmap
     drawn from the pooled-across-replicates ``*_pooled_contact_map.npy``
     when one is present.
     File: ``{condition}_contact_map_merged_reps.png``

  3. **All-conditions overview** — one subplot per condition (using the
     pooled map when available, else replicate 0), showing every condition
     side-by-side on a shared scale.
     File: ``all_conditions_contact_maps.png``

The script discovers what's available by globbing ``results/analysis/``;
it does not need to know the list of conditions in advance.

Usage::

    python scripts/plot_contact_map_panels.py
    python scripts/plot_contact_map_panels.py --analysis-dir results/analysis \\
                                              --output-dir results/figures
"""

from __future__ import annotations

import argparse
import glob
import logging
import math
import os
import re
import sys
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# Filename patterns produced by run_analysis_all.py (Tasks 5+6 of the
# 2026-04-27 flat-output refactor):
#   {condition}_{n_blocks}blk_rep{rep}_contact_map.npy
#   {condition}_{total_blocks}blk_pooled_contact_map.npy
PER_REP_RE = re.compile(r"^(?P<cond>.+?)_(?P<blocks>\d+)blk_rep(?P<rep>\d+)_contact_map\.npy$")
POOLED_RE = re.compile(r"^(?P<cond>.+?)_(?P<blocks>\d+)blk_pooled_contact_map\.npy$")


def discover_maps(analysis_dir: str) -> tuple[dict, dict]:
    """Return (per_rep, pooled) dicts.

    per_rep[condition] -> list of (rep:int, path:str) sorted by rep.
    pooled[condition]  -> (n_blocks:int, path:str)  (latest by n_blocks).
    """
    per_rep: dict[str, list[tuple[int, str]]] = defaultdict(list)
    pooled: dict[str, tuple[int, str]] = {}

    for path in sorted(glob.glob(os.path.join(analysis_dir, "*_contact_map.npy"))):
        fn = os.path.basename(path)
        m = POOLED_RE.match(fn)
        if m:
            cond = m.group("cond")
            blocks = int(m.group("blocks"))
            # If multiple pooled files exist for the same condition keep
            # the one with the largest block count (latest run).
            if cond not in pooled or blocks > pooled[cond][0]:
                pooled[cond] = (blocks, path)
            continue
        m = PER_REP_RE.match(fn)
        if m:
            cond = m.group("cond")
            rep = int(m.group("rep"))
            per_rep[cond].append((rep, path))

    for cond in per_rep:
        per_rep[cond].sort(key=lambda x: x[0])

    return per_rep, pooled


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _import_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    return plt, LogNorm


def _plot_one(ax, matrix, title, vmin, vmax, log_scale: bool):
    plt, LogNorm = _import_mpl()
    import numpy as np
    data = matrix.astype(np.float64, copy=False)
    if log_scale:
        # Floor zeros to a small positive so LogNorm doesn't choke.
        floor = max(np.nanmin(data[data > 0]) if np.any(data > 0) else 1e-6, 1e-12)
        plot_data = np.where(data > 0, data, floor)
        norm = LogNorm(vmin=max(vmin, floor), vmax=vmax)
    else:
        plot_data = data
        norm = None
    im = ax.imshow(plot_data, cmap="Reds", norm=norm,
                   origin="upper", interpolation="nearest")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    return im


def _shared_scale(matrices, percentile: float = 98.0):
    import numpy as np
    flat = np.concatenate([m[m > 0].ravel() for m in matrices if (m > 0).any()])
    if flat.size == 0:
        return 1e-6, 1.0
    vmin = float(np.percentile(flat, 1.0))
    vmax = float(np.percentile(flat, percentile))
    if vmax <= vmin:
        vmax = vmin * 10 if vmin > 0 else 1.0
    return vmin, vmax


def _grid_dims(n: int) -> tuple[int, int]:
    """Pick a (rows, cols) layout that's roughly square for n subplots."""
    cols = int(math.ceil(math.sqrt(n)))
    rows = int(math.ceil(n / cols))
    return rows, cols


# ---------------------------------------------------------------------------
# Panel builders
# ---------------------------------------------------------------------------

def plot_per_rep_panel(condition: str, rep_paths: list[tuple[int, str]],
                       output_path: str, log_scale: bool = True) -> None:
    """One subplot per replicate of a single condition."""
    plt, _ = _import_mpl()
    import numpy as np

    if not rep_paths:
        logger.warning(f"  [{condition}] no replicate maps; skipping per-rep panel")
        return

    matrices = [np.load(p) for _, p in rep_paths]
    vmin, vmax = _shared_scale(matrices)

    n = len(matrices)
    rows, cols = _grid_dims(n)
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows),
                             squeeze=False)
    last_im = None
    for i, ((rep, path), m) in enumerate(zip(rep_paths, matrices)):
        ax = axes[i // cols][i % cols]
        last_im = _plot_one(ax, m, f"rep{rep}", vmin, vmax, log_scale)
    # Blank any leftover axes.
    for j in range(n, rows * cols):
        axes[j // cols][j % cols].axis("off")

    fig.suptitle(f"{condition}: per-replicate contact maps", fontsize=13)
    if last_im is not None:
        fig.colorbar(last_im, ax=axes.ravel().tolist(),
                     fraction=0.025, pad=0.02, label="contact freq")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  wrote {output_path}")


def plot_pooled_panel(condition: str, pooled_path: str,
                      output_path: str, log_scale: bool = True) -> None:
    """Single contact-map heatmap from the pooled-across-reps file."""
    plt, _ = _import_mpl()
    import numpy as np

    matrix = np.load(pooled_path)
    vmin, vmax = _shared_scale([matrix])

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    im = _plot_one(ax, matrix, f"{condition} (pooled)", vmin, vmax, log_scale)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="contact freq")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  wrote {output_path}")


def plot_all_conditions_panel(per_cond_map: dict[str, "np.ndarray"],
                              output_path: str, log_scale: bool = True) -> None:
    """One subplot per condition, side by side, on a shared colour scale."""
    plt, _ = _import_mpl()

    if not per_cond_map:
        logger.warning("  [all_conditions] no maps; skipping overview panel")
        return

    items = sorted(per_cond_map.items())
    matrices = [m for _, m in items]
    vmin, vmax = _shared_scale(matrices)

    n = len(items)
    rows, cols = _grid_dims(n)
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows),
                             squeeze=False)
    last_im = None
    for i, ((cond, m)) in enumerate(items):
        ax = axes[i // cols][i % cols]
        last_im = _plot_one(ax, m, cond, vmin, vmax, log_scale)
    for j in range(n, rows * cols):
        axes[j // cols][j % cols].axis("off")

    fig.suptitle("All conditions: contact maps", fontsize=13)
    if last_im is not None:
        fig.colorbar(last_im, ax=axes.ravel().tolist(),
                     fraction=0.025, pad=0.02, label="contact freq")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  wrote {output_path}")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--analysis-dir", default="results/analysis",
                        help="Folder holding *_contact_map.npy files. "
                             "Default: results/analysis.")
    parser.add_argument("--output-dir", default="results/figures",
                        help="Destination for the figure panels. "
                             "Default: results/figures.")
    parser.add_argument("--no-log-scale", action="store_true",
                        help="Plot raw contact frequencies instead of log "
                             "(default: log scale, better for Hi-C contact maps).")
    args = parser.parse_args()

    if not os.path.isdir(args.analysis_dir):
        logger.error(f"--analysis-dir does not exist: {args.analysis_dir}")
        return 1
    os.makedirs(args.output_dir, exist_ok=True)

    import numpy as np  # noqa: F401  (used inside helpers via lazy imports)

    per_rep, pooled = discover_maps(args.analysis_dir)
    if not per_rep and not pooled:
        logger.error(
            f"No contact-map files found under {args.analysis_dir}. "
            f"Expected files matching '<condition>_<blocks>blk_rep<N>_contact_map.npy' "
            f"or '<condition>_<blocks>blk_pooled_contact_map.npy'."
        )
        return 1

    log_scale = not args.no_log_scale
    conditions = sorted(set(per_rep) | set(pooled))
    logger.info(f"Discovered {len(conditions)} condition(s): {conditions}")

    # 1. Per-condition per-rep grids.
    for cond in conditions:
        if cond not in per_rep:
            continue
        plot_per_rep_panel(
            cond, per_rep[cond],
            os.path.join(args.output_dir, f"{cond}_contact_maps_per_rep.png"),
            log_scale=log_scale,
        )

    # 2. Per-condition pooled (merged-reps).
    for cond in conditions:
        if cond not in pooled:
            logger.info(f"  [{cond}] no pooled map; skipping merged-reps panel")
            continue
        _, path = pooled[cond]
        plot_pooled_panel(
            cond, path,
            os.path.join(args.output_dir, f"{cond}_contact_map_merged_reps.png"),
            log_scale=log_scale,
        )

    # 3. All-conditions overview.
    import numpy as np
    overview: dict[str, np.ndarray] = {}
    for cond in conditions:
        if cond in pooled:
            overview[cond] = np.load(pooled[cond][1])
        elif cond in per_rep:
            # Fall back to rep0 (the first replicate available).
            overview[cond] = np.load(per_rep[cond][0][1])
    plot_all_conditions_panel(
        overview, os.path.join(args.output_dir, "all_conditions_contact_maps.png"),
        log_scale=log_scale,
    )

    logger.info(f"All panels written to {args.output_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
