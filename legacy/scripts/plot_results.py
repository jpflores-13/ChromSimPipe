#!/usr/bin/env python
"""
Visualization of simulation results and comparison with experimental Hi-C.

Generates publication-quality figures:
  1. Contact map heatmaps (simulated vs experimental)
  2. P(s) curves overlaid
  3. Insulation score profiles
  4. Parameter sweep summary (which neuron parameters best match data)

Usage:
    python plot_results.py --results-dir ../results --data-dir ../data --output ../results/figures
"""

import os
import sys
import argparse
import json
import logging

if os.environ.get("CONDA_DEFAULT_ENV", "") != "cohesin_sim":
    sys.stderr.write(
        f"[env] active conda env is '{os.environ.get('CONDA_DEFAULT_ENV') or 'none'}', "
        "expected 'cohesin_sim'. run: conda activate cohesin_sim\n"
    )
    sys.exit(1)

import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import (
    N_MONOMERS, RESOLUTION, REGION_START, REGION_END,
    SOX2_START_MONOMER, SOX2_END_MONOMER, ALL_PARAM_SETS,
)


def plot_contact_map(
    matrix: np.ndarray,
    title: str,
    output_path: str,
    vmax_percentile: float = 98,
    region_start: int = REGION_START,
    resolution: int = RESOLUTION,
    sox2_pos: int = SOX2_START_MONOMER,
    log_scale: bool = True,
):
    """Plot a single contact map heatmap."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))

    # Log transform for better visualization
    plot_data = matrix.copy()
    if log_scale:
        plot_data = np.log2(plot_data + 1)

    vmax = np.percentile(plot_data[plot_data > 0], vmax_percentile) if (plot_data > 0).any() else 1
    vmin = 0

    # Genomic coordinates for axis labels
    n = matrix.shape[0]
    extent_mb = [
        region_start / 1e6,
        (region_start + n * resolution) / 1e6,
        (region_start + n * resolution) / 1e6,
        region_start / 1e6,
    ]

    im = ax.imshow(plot_data, cmap="fall" if "fall" in plt.colormaps() else "YlOrRd",
                   vmin=vmin, vmax=vmax, extent=extent_mb, aspect="equal",
                   interpolation="none")

    # Mark Sox2 position
    sox2_mb = (region_start + sox2_pos * resolution) / 1e6
    ax.axhline(sox2_mb, color="blue", linewidth=0.5, alpha=0.5, linestyle="--")
    ax.axvline(sox2_mb, color="blue", linewidth=0.5, alpha=0.5, linestyle="--")
    ax.text(sox2_mb, extent_mb[0] + 0.05, "Sox2", fontsize=8, color="blue", ha="center")

    ax.set_xlabel(f"Genomic position (Mb) — chr3", fontsize=12)
    ax.set_ylabel(f"Genomic position (Mb) — chr3", fontsize=12)
    ax.set_title(title, fontsize=14)

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log2(contact frequency + 1)" if log_scale else "Contact frequency")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_contact_map_comparison(
    maps: dict,
    output_path: str,
    titles: dict = None,
):
    """
    Plot multiple contact maps side by side for comparison.

    Parameters
    ----------
    maps : dict
        {name: contact_matrix} pairs
    titles : dict, optional
        {name: display_title} pairs
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_maps = len(maps)
    fig, axes = plt.subplots(1, n_maps, figsize=(6 * n_maps, 5.5))
    if n_maps == 1:
        axes = [axes]

    for ax, (name, matrix) in zip(axes, maps.items()):
        plot_data = np.log2(matrix + 1)
        vmax = np.percentile(plot_data[plot_data > 0], 98) if (plot_data > 0).any() else 1

        n = matrix.shape[0]
        extent_mb = [
            REGION_START / 1e6,
            (REGION_START + n * RESOLUTION) / 1e6,
            (REGION_START + n * RESOLUTION) / 1e6,
            REGION_START / 1e6,
        ]

        im = ax.imshow(plot_data, cmap="YlOrRd", vmin=0, vmax=vmax,
                       extent=extent_mb, aspect="equal", interpolation="none")

        title = titles.get(name, name) if titles else name
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Position (Mb)")
        if ax == axes[0]:
            ax.set_ylabel("Position (Mb)")

        # Mark Sox2
        sox2_mb = (REGION_START + SOX2_START_MONOMER * RESOLUTION) / 1e6
        ax.axhline(sox2_mb, color="blue", lw=0.5, alpha=0.4, ls="--")
        ax.axvline(sox2_mb, color="blue", lw=0.5, alpha=0.4, ls="--")

    plt.suptitle("Contact Maps — Sox2 Locus (chr3:34-36 Mb, mm10)", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_ps_curves(
    ps_data: dict,
    output_path: str,
    exp_ps: dict = None,
):
    """
    Plot P(s) curves for multiple conditions.

    Parameters
    ----------
    ps_data : dict
        {name: (distances, ps)} from simulations
    exp_ps : dict, optional
        {name: (distances, ps)} from experimental Hi-C
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    colors_sim = plt.cm.tab10(np.linspace(0, 1, len(ps_data)))
    colors_exp = ["black", "gray"]

    # Plot experimental first
    if exp_ps:
        for idx, (name, (dist, ps)) in enumerate(exp_ps.items()):
            mask = (dist > 0) & (ps > 0)
            ax.loglog(dist[mask] / 1e3, ps[mask],
                      color=colors_exp[idx % len(colors_exp)],
                      linewidth=2.5, label=f"Exp: {name}", linestyle="-")

    # Plot simulated
    for idx, (name, (dist, ps)) in enumerate(ps_data.items()):
        mask = (dist > 0) & (ps > 0)
        ax.loglog(dist[mask] / 1e3, ps[mask],
                  color=colors_sim[idx], linewidth=1.5,
                  label=f"Sim: {name}", linestyle="--")

    ax.set_xlabel("Genomic distance (kb)", fontsize=12)
    ax.set_ylabel("Contact probability P(s)", fontsize=12)
    ax.set_title("Contact Probability Decay — Sox2 Locus", fontsize=13)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1, 2000)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_insulation_profiles(
    insulation_data: dict,
    output_path: str,
    exp_insulation: dict = None,
):
    """
    Plot insulation score profiles for multiple conditions.

    Parameters
    ----------
    insulation_data : dict
        {name: insulation_scores} from simulations
    exp_insulation : dict, optional
        {name: insulation_scores} from experimental Hi-C
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))

    positions_mb = (REGION_START + np.arange(N_MONOMERS) * RESOLUTION) / 1e6

    # Experimental
    if exp_insulation:
        for name, ins in exp_insulation.items():
            n = min(len(ins), len(positions_mb))
            ax.plot(positions_mb[:n], ins[:n], color="black", linewidth=2,
                    label=f"Exp: {name}", alpha=0.8)

    # Simulated
    colors = plt.cm.tab10(np.linspace(0, 1, len(insulation_data)))
    for idx, (name, ins) in enumerate(insulation_data.items()):
        n = min(len(ins), len(positions_mb))
        ax.plot(positions_mb[:n], ins[:n], color=colors[idx],
                linewidth=1.2, label=f"Sim: {name}", alpha=0.7)

    # Mark Sox2
    sox2_mb = (REGION_START + SOX2_START_MONOMER * RESOLUTION) / 1e6
    ax.axvline(sox2_mb, color="blue", lw=1, alpha=0.5, ls="--", label="Sox2")

    ax.set_xlabel("Genomic position (Mb) — chr3", fontsize=12)
    ax.set_ylabel("Insulation score (log2)", fontsize=12)
    ax.set_title("Insulation Score Profile — Sox2 Locus", fontsize=13)
    ax.legend(fontsize=7, loc="lower right", ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


def plot_parameter_sweep_summary(
    metrics_per_condition: dict,
    output_path: str,
):
    """
    Bar chart comparing how well each parameter set matches experimental data.

    Parameters
    ----------
    metrics_per_condition : dict
        {condition_name: metrics_dict} from compare_contact_maps()
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conditions = list(metrics_per_condition.keys())
    metric_names = ["overall_pearson", "stratum_adjusted_corr", "ps_curve_corr", "insulation_corr"]
    metric_labels = ["Overall Pearson", "Stratum-adj. Corr", "P(s) Corr", "Insulation Corr"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()

    x = np.arange(len(conditions))
    width = 0.6

    for idx, (mname, mlabel) in enumerate(zip(metric_names, metric_labels)):
        values = [metrics_per_condition[c].get(mname, 0) for c in conditions]

        bars = axes[idx].bar(x, values, width, color=plt.cm.Set2(np.linspace(0, 1, len(conditions))))
        axes[idx].set_ylabel(mlabel)
        axes[idx].set_xticks(x)
        axes[idx].set_xticklabels([c.replace("_", "\n") for c in conditions],
                                   fontsize=7, rotation=45, ha="right")
        axes[idx].set_ylim(0, 1)
        axes[idx].axhline(0.5, color="gray", ls="--", alpha=0.3)

        # Annotate best
        best_idx = np.argmax(values)
        axes[idx].get_children()[best_idx].set_edgecolor("red")
        axes[idx].get_children()[best_idx].set_linewidth(2)

    plt.suptitle("Parameter Sweep: Similarity to Experimental Hi-C\n(Sox2 locus, Cortical Neurons)",
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Plot simulation results")
    parser.add_argument("--results-dir", type=str, required=True,
                        help="Directory containing simulation results")
    parser.add_argument("--data-dir", type=str, default=None,
                        help="Directory containing experimental Hi-C .npy files")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory for figures")

    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(args.results_dir, "figures")
    os.makedirs(args.output, exist_ok=True)

    # --- Collect all simulation results ---
    sim_maps = {}
    sim_ps = {}
    sim_insulation = {}
    metrics_vs_cn = {}

    for param_set in ALL_PARAM_SETS:
        name = param_set["name"]
        for rep in range(3):  # check replicates
            analysis_dir = os.path.join(args.results_dir, f"{name}_rep{rep}", "analysis")
            if not os.path.exists(analysis_dir):
                continue

            map_path = os.path.join(analysis_dir, "sim_contact_map.npy")
            if os.path.exists(map_path):
                sim_maps[f"{name}_r{rep}"] = np.load(map_path)

            ps_path = os.path.join(analysis_dir, "sim_ps_curve.npz")
            if os.path.exists(ps_path):
                data = np.load(ps_path)
                sim_ps[f"{name}_r{rep}"] = (data["distances"], data["ps"])

            ins_path = os.path.join(analysis_dir, "sim_insulation.npy")
            if os.path.exists(ins_path):
                sim_insulation[f"{name}_r{rep}"] = np.load(ins_path)

            metrics_path = os.path.join(analysis_dir, "comparison_metrics.json")
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    metrics_vs_cn[f"{name}_r{rep}"] = json.load(f)

    if not sim_maps:
        logger.warning("No simulation results found. Run simulations first.")
        logger.info("Expected structure: results/<name>_rep<N>/analysis/sim_contact_map.npy")
        return

    # --- Load experimental data ---
    exp_maps = {}
    exp_ps = {}
    exp_insulation = {}

    if args.data_dir:
        for cell_type in ["mESC", "CN"]:
            map_path = os.path.join(args.data_dir, f"hic_{cell_type}_Sox2.npy")
            if os.path.exists(map_path):
                exp_maps[cell_type] = np.load(map_path)

            ins_path = os.path.join(args.data_dir, f"insulation_{cell_type}.npy")
            if os.path.exists(ins_path):
                exp_insulation[cell_type] = np.load(ins_path)

    # --- Generate plots ---
    # Individual contact maps
    for name, cmap in sim_maps.items():
        plot_contact_map(cmap, f"Simulated: {name}",
                        os.path.join(args.output, f"cmap_{name}.png"))

    for name, cmap in exp_maps.items():
        plot_contact_map(cmap, f"Experimental: {name}",
                        os.path.join(args.output, f"cmap_exp_{name}.png"))

    # Comparison panel
    if sim_maps:
        # Pick one representative per condition (rep0)
        repr_maps = {}
        for ps in ALL_PARAM_SETS:
            key = f"{ps['name']}_r0"
            if key in sim_maps:
                repr_maps[ps['name']] = sim_maps[key]
        repr_maps.update(exp_maps)
        if repr_maps:
            plot_contact_map_comparison(repr_maps,
                                       os.path.join(args.output, "cmap_comparison.png"))

    # P(s) curves
    if sim_ps:
        plot_ps_curves(sim_ps, os.path.join(args.output, "ps_curves.png"),
                       exp_ps=exp_ps if exp_ps else None)

    # Insulation profiles
    if sim_insulation:
        plot_insulation_profiles(sim_insulation,
                                os.path.join(args.output, "insulation_profiles.png"),
                                exp_insulation=exp_insulation if exp_insulation else None)

    # Parameter sweep summary
    if metrics_vs_cn:
        plot_parameter_sweep_summary(metrics_vs_cn,
                                    os.path.join(args.output, "parameter_sweep.png"))

    logger.info(f"All figures saved to: {args.output}")


if __name__ == "__main__":
    main()
