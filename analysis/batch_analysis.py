#!/usr/bin/env python
"""
Batch analysis: compute contact maps for all completed conditions and
generate a multi-panel comparison figure.

Scans results/polychrom_3d/ for completed simulation directories,
runs contact map analysis on each, saves individual PNGs per condition,
and produces a single multi-panel figure for quick comparison.

Usage:
    python analysis/batch_analysis.py
    python analysis/batch_analysis.py --results-dir results/polychrom_3d --n-jobs 4
    python analysis/batch_analysis.py --skip-existing   # only analyze new results
"""

import os
import sys
import glob
import argparse
import logging
import json

if os.environ.get("CONDA_DEFAULT_ENV", "") != "cohesin_sim":
    sys.stderr.write(
        f"[env] active conda env is '{os.environ.get('CONDA_DEFAULT_ENV') or 'none'}', "
        "expected 'cohesin_sim'. run: conda activate cohesin_sim\n"
    )
    sys.exit(1)

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import (
    N_MONOMERS, RESOLUTION, REGION_START, REGION_END,
    SOX2_START_MONOMER, SOX2_END_MONOMER,
    SIMULATION_CONDITIONS,
)
from analysis.contact_maps import (
    load_conformations_h5,
    compute_contact_map_from_conformations,
    extract_tiles_and_average,
    compute_ps_curve,
    compute_insulation_score,
    load_lef_contact_map,
)


def find_completed_dirs(results_dir):
    """
    Scan results directory for completed simulation runs.
    Returns list of (condition_label, dir_path) sorted by condition name.
    """
    completed = []
    if not os.path.isdir(results_dir):
        return completed

    for entry in sorted(os.listdir(results_dir)):
        full = os.path.join(results_dir, entry)
        if not os.path.isdir(full):
            continue
        # Skip shard directories (they have _shardN suffix)
        if "_shard" in entry:
            continue
        # Check it has actual simulation output
        has_blocks = bool(glob.glob(os.path.join(full, "blocks_*.h5")))
        has_conf = os.path.exists(os.path.join(full, "conformations.h5"))
        has_lef = os.path.exists(os.path.join(full, "lef_contact_map.npy"))
        if has_blocks or has_conf or has_lef:
            completed.append((entry, full))

    return completed


def run_analysis_for_dir(sim_dir, contact_radius=3.0, n_jobs=4, skip_existing=False):
    """
    Run contact map analysis for a single simulation directory.
    Returns the contact map numpy array or None on failure.
    """
    analysis_dir = os.path.join(sim_dir, "analysis")
    map_path = os.path.join(analysis_dir, "sim_contact_map.npy")

    if skip_existing and os.path.exists(map_path):
        logger.info(f"  Loading existing contact map")
        return np.load(map_path)

    os.makedirs(analysis_dir, exist_ok=True)

    # Load conformations
    conformations = None
    try:
        conformations = load_conformations_h5(sim_dir)
    except FileNotFoundError:
        pass

    if conformations is not None:
        logger.info(f"  Computing contact map from {len(conformations)} frames "
                    f"({n_jobs} workers)...")
        # Use tile extraction (handles both tiled and single-locus conformations)
        sim_map = extract_tiles_and_average(
            conformations, contact_radius, n_jobs=n_jobs)
    else:
        sim_map = load_lef_contact_map(sim_dir)
        if sim_map is None:
            logger.warning(f"  No data found, skipping")
            return None
        logger.info(f"  Using LEF bridging contact map")

    # Save
    np.save(map_path, sim_map)
    distances, ps = compute_ps_curve(sim_map, RESOLUTION)
    np.savez(os.path.join(analysis_dir, "sim_ps_curve.npz"),
             distances=distances, ps=ps)
    insulation = compute_insulation_score(sim_map)
    np.save(os.path.join(analysis_dir, "sim_insulation.npy"), insulation)

    return sim_map


def plot_single_contact_map(contact_map, title, output_path):
    """Save a single contact map as PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 6))

    # Genomic coordinates in Mb
    extent_mb = [
        REGION_START / 1e6, REGION_END / 1e6,
        REGION_END / 1e6, REGION_START / 1e6,
    ]

    im = ax.imshow(
        np.log2(contact_map + 1e-6),
        cmap="hot_r", vmin=-10, vmax=0,
        origin="upper", extent=extent_mb,
    )

    # Mark Sox2 locus
    sox2_start_mb = (REGION_START + SOX2_START_MONOMER * RESOLUTION) / 1e6
    sox2_end_mb = (REGION_START + SOX2_END_MONOMER * RESOLUTION) / 1e6
    for pos in [sox2_start_mb, sox2_end_mb]:
        ax.axhline(pos, color="cyan", lw=0.5, alpha=0.5)
        ax.axvline(pos, color="cyan", lw=0.5, alpha=0.5)

    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Position (Mb)")
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, label="log2(contact freq)", shrink=0.8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_multipanel_figure(results, output_path):
    """
    Create a multi-panel figure with one contact map per condition.

    Parameters
    ----------
    results : list of (label, contact_map)
    output_path : str
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(results)
    if n == 0:
        logger.warning("No results to plot")
        return

    # Layout: aim for roughly square grid
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    if n == 1:
        axes = np.array([axes])
    axes = np.atleast_2d(axes)

    extent_mb = [
        REGION_START / 1e6, REGION_END / 1e6,
        REGION_END / 1e6, REGION_START / 1e6,
    ]

    sox2_start_mb = (REGION_START + SOX2_START_MONOMER * RESOLUTION) / 1e6
    sox2_end_mb = (REGION_START + SOX2_END_MONOMER * RESOLUTION) / 1e6

    for idx, (label, cmap) in enumerate(results):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        im = ax.imshow(
            np.log2(cmap + 1e-6),
            cmap="hot_r", vmin=-10, vmax=0,
            origin="upper", extent=extent_mb,
        )

        # Mark Sox2
        for pos in [sox2_start_mb, sox2_end_mb]:
            ax.axhline(pos, color="cyan", lw=0.3, alpha=0.4)
            ax.axvline(pos, color="cyan", lw=0.3, alpha=0.4)

        # Clean up label for title
        title = label.replace("_neuron_ctcf", "\n(neuron CTCF)")
        title = title.replace("_ctcf-mESC", "\n(mESC CTCF)")
        title = title.replace("_ctcf-neuron", "\n(neuron CTCF)")
        ax.set_title(title, fontsize=8, fontweight="bold")
        ax.set_xlabel("Mb", fontsize=7)
        ax.set_ylabel("Mb", fontsize=7)
        ax.tick_params(labelsize=6)

    # Turn off empty panels
    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].axis("off")

    # Shared colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="log2(contact freq)")

    fig.suptitle("Simulated contact maps — Sox2 locus (chr3:34–36 Mb, mm10)",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Multi-panel figure saved: {output_path}")


def plot_ps_comparison(results_dirs, output_path):
    """
    Overlay P(s) curves from all conditions on one plot.

    Parameters
    ----------
    results_dirs : list of (label, sim_dir)
    output_path : str
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))

    for label, sim_dir in results_dirs:
        ps_path = os.path.join(sim_dir, "analysis", "sim_ps_curve.npz")
        if not os.path.exists(ps_path):
            continue
        data = np.load(ps_path)
        d, p = data["distances"], data["ps"]
        mask = (d > 0) & (p > 0)
        short_label = label.split("_ctcf-")[0] if "_ctcf-" in label else label
        ax.plot(d[mask] / 1e3, p[mask], label=short_label, lw=1.2, alpha=0.8)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genomic distance (kb)")
    ax.set_ylabel("Contact probability P(s)")
    ax.set_title("P(s) curves — all conditions")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logger.info(f"P(s) comparison saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Batch analysis of all completed conditions")
    parser.add_argument("--results-dir", type=str, default="results/polychrom_3d",
                        help="Base results directory")
    parser.add_argument("--output-dir", type=str, default="results/figures",
                        help="Output directory for PNGs")
    parser.add_argument("--n-jobs", type=int, default=4,
                        help="Parallel workers for contact map computation")
    parser.add_argument("--contact-radius", type=float, default=3.0,
                        help="Contact detection radius")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip conditions that already have contact maps")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Find completed simulation directories ---
    completed = find_completed_dirs(args.results_dir)
    if not completed:
        logger.error(f"No completed simulations found in {args.results_dir}")
        sys.exit(1)

    logger.info(f"Found {len(completed)} completed simulation(s):")
    for label, path in completed:
        logger.info(f"  {label}")

    # --- Analyze each condition ---
    panel_data = []  # (label, contact_map) for multipanel figure

    for i, (label, sim_dir) in enumerate(completed):
        logger.info(f"\n[{i+1}/{len(completed)}] Analyzing: {label}")

        cmap = run_analysis_for_dir(
            sim_dir,
            contact_radius=args.contact_radius,
            n_jobs=args.n_jobs,
            skip_existing=args.skip_existing,
        )

        if cmap is not None:
            # Save individual PNG
            png_path = os.path.join(args.output_dir, f"{label}.png")
            plot_single_contact_map(cmap, label, png_path)
            logger.info(f"  Saved: {png_path}")

            panel_data.append((label, cmap))

    # --- Multi-panel comparison figure ---
    if panel_data:
        multipanel_path = os.path.join(args.output_dir, "all_conditions_comparison.png")
        plot_multipanel_figure(panel_data, multipanel_path)

        # P(s) overlay
        ps_path = os.path.join(args.output_dir, "ps_curves_comparison.png")
        plot_ps_comparison(completed, ps_path)

    logger.info(f"\nDone! All figures saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
