#!/usr/bin/env python
"""
Batch runner: launch all parameter conditions and replicates.

Usage:
    # Run all conditions, 3 replicates each (LEF-only, no GPU needed)
    python run_batch.py --engine lef_only

    # Run all conditions with polychrom (requires GPU)
    python run_batch.py --engine polychrom --gpu 0

    # Run specific conditions
    python run_batch.py --conditions mESC CN_long_residency --engine lef_only

    # Then analyze all
    python run_batch.py --analyze-only --exp-es ../data/hic_mESC_Sox2.npy --exp-cn ../data/hic_CN_Sox2.npy
"""

import os
import sys
import argparse
import logging
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import ALL_PARAM_SETS, SIMULATION_CONDITIONS, SIM_RUN
from scripts.run_simulation import (
    run_simulation_polychrom, run_simulation_openmm_standalone,
    get_param_set, resolve_condition,
)
from analysis.contact_maps import (
    load_conformations_h5, load_lef_contact_map,
    compute_contact_map_from_conformations, compute_ps_curve,
    compute_insulation_score, compare_contact_maps,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def run_all_simulations(conditions, n_replicates, engine, output_dir, gpu=0):
    """
    Run simulations for all specified conditions and replicates.

    Parameters
    ----------
    conditions : list of dict
        Each dict is a SIMULATION_CONDITIONS entry with keys:
        name, params, ctcf_type, compare_hic.
    n_replicates : int
    engine : str
    output_dir : str
    gpu : int
    """
    results = {}

    for cond in conditions:
        cond_name = cond["name"]
        params = cond["params"]
        ctcf_type = cond["ctcf_type"]

        for rep in range(n_replicates):
            logger.info(f"\n{'='*60}")
            logger.info(f"RUNNING: {cond_name} (params={params['name']}, "
                        f"CTCF={ctcf_type}), replicate {rep}")
            logger.info(f"{'='*60}")

            run_kwargs = dict(params=params, replicate=rep,
                              output_dir=output_dir, gpu=gpu,
                              ctcf_type=ctcf_type)

            if engine == "polychrom":
                run_dir = run_simulation_polychrom(**run_kwargs)
            elif engine == "openmm":
                run_dir = run_simulation_openmm_standalone(**run_kwargs)
            else:
                run_dir = run_simulation_polychrom(**run_kwargs)

            results[f"{cond_name}_rep{rep}"] = run_dir

    return results


def analyze_all_results(output_dir, exp_es_path=None, exp_cn_path=None):
    """Run analysis on all simulation results."""
    from configs.parameters import RESOLUTION, SIM_RUN

    # Load experimental data if available
    exp_es = np.load(exp_es_path) if exp_es_path and os.path.exists(exp_es_path) else None
    exp_cn = np.load(exp_cn_path) if exp_cn_path and os.path.exists(exp_cn_path) else None

    all_metrics = {}

    for cond in SIMULATION_CONDITIONS:
        cond_name = cond["name"]
        params = cond["params"]
        ctcf_type = cond["ctcf_type"]
        compare_hic = cond.get("compare_hic")
        name = params["name"]

        for rep in range(SIM_RUN["n_replicates"]):
            run_dir = os.path.join(output_dir,
                                   f"{name}_ctcf-{ctcf_type}_rep{rep}")
            if not os.path.exists(run_dir):
                continue

            analysis_dir = os.path.join(run_dir, "analysis")
            os.makedirs(analysis_dir, exist_ok=True)

            logger.info(f"\nAnalyzing: {cond_name}_rep{rep}")

            # Load simulation data
            try:
                conformations = load_conformations_h5(run_dir)
                if conformations:
                    sim_map = compute_contact_map_from_conformations(
                        conformations, SIM_RUN["contact_radius"])
                else:
                    sim_map = load_lef_contact_map(run_dir)
            except FileNotFoundError:
                sim_map = load_lef_contact_map(run_dir)

            if sim_map is None:
                logger.warning(f"  No data found, skipping")
                continue

            # Save contact map
            np.save(os.path.join(analysis_dir, "sim_contact_map.npy"), sim_map)

            # P(s) curve
            distances, ps = compute_ps_curve(sim_map, RESOLUTION)
            np.savez(os.path.join(analysis_dir, "sim_ps_curve.npz"),
                     distances=distances, ps=ps)

            # Insulation score
            insulation = compute_insulation_score(sim_map)
            np.save(os.path.join(analysis_dir, "sim_insulation.npy"), insulation)

            # Compare with the Hi-C dataset this condition should match
            key = f"{cond_name}_rep{rep}"

            # Always compare against both Hi-C datasets (when available)
            # but flag which one is the "expected" match
            if exp_es is not None:
                metrics_es = compare_contact_maps(sim_map, exp_es)
                metrics_es["is_target"] = (compare_hic == "mESC")
                with open(os.path.join(analysis_dir, "metrics_vs_mESC.json"), "w") as f:
                    json.dump(metrics_es, f, indent=2)
                all_metrics[f"{key}_vs_mESC"] = metrics_es
                marker = " ← target" if compare_hic == "mESC" else ""
                logger.info(f"  vs mESC: Pearson={metrics_es['overall_pearson']:.3f}, "
                           f"SCC={metrics_es['stratum_adjusted_corr']:.3f}{marker}")

            if exp_cn is not None:
                metrics_cn = compare_contact_maps(sim_map, exp_cn)
                metrics_cn["is_target"] = (compare_hic == "neuron")
                with open(os.path.join(analysis_dir, "metrics_vs_CN.json"), "w") as f:
                    json.dump(metrics_cn, f, indent=2)
                all_metrics[f"{key}_vs_CN"] = metrics_cn
                marker = " ← target" if compare_hic == "neuron" else ""
                logger.info(f"  vs CN:   Pearson={metrics_cn['overall_pearson']:.3f}, "
                           f"SCC={metrics_cn['stratum_adjusted_corr']:.3f}{marker}")

    # Save summary
    summary_path = os.path.join(output_dir, "analysis_summary.json")
    with open(summary_path, "w") as f:
        json.dump(all_metrics, f, indent=2)
    logger.info(f"\nSummary saved to: {summary_path}")

    # Print summary table
    if all_metrics:
        logger.info("\n" + "="*100)
        logger.info("SUMMARY: Similarity to Experimental Hi-C")
        logger.info("="*100)
        logger.info(f"{'Condition':<42} {'CTCF':>6} {'Target':>7}  "
                    f"{'vs mESC (SCC)':<15} {'vs CN (SCC)':<15}")
        logger.info("-"*95)

        for cond in SIMULATION_CONDITIONS:
            cond_name = cond["name"]
            ctcf_type = cond["ctcf_type"]
            target = cond.get("compare_hic") or "—"

            scc_es_vals = []
            scc_cn_vals = []
            for rep in range(SIM_RUN["n_replicates"]):
                key_es = f"{cond_name}_rep{rep}_vs_mESC"
                key_cn = f"{cond_name}_rep{rep}_vs_CN"
                if key_es in all_metrics:
                    scc_es_vals.append(all_metrics[key_es]["stratum_adjusted_corr"])
                if key_cn in all_metrics:
                    scc_cn_vals.append(all_metrics[key_cn]["stratum_adjusted_corr"])

            scc_es = f"{np.mean(scc_es_vals):.3f}±{np.std(scc_es_vals):.3f}" if scc_es_vals else "N/A"
            scc_cn = f"{np.mean(scc_cn_vals):.3f}±{np.std(scc_cn_vals):.3f}" if scc_cn_vals else "N/A"
            logger.info(f"{cond_name:<42} {ctcf_type:>6} {target:>7}  "
                        f"{scc_es:<15} {scc_cn:<15}")

    return all_metrics


def main():
    parser = argparse.ArgumentParser(description="Batch simulation runner")
    parser.add_argument("--conditions", nargs="+", default=None,
                        help="Condition names to run from SIMULATION_CONDITIONS "
                             "(e.g., mESC_ctrl CN_long_residency_neuron_ctcf). "
                             "Default: all conditions.")
    parser.add_argument("--n-replicates", type=int, default=SIM_RUN["n_replicates"],
                        help="Number of replicates per condition")
    parser.add_argument("--engine", type=str, default="lef_only",
                        choices=["polychrom", "openmm", "lef_only"],
                        help="Simulation engine")
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--analyze-only", action="store_true",
                        help="Skip simulations, only analyze existing results")
    parser.add_argument("--exp-es", type=str, default=None,
                        help="Path to mESC experimental Hi-C .npy")
    parser.add_argument("--exp-cn", type=str, default=None,
                        help="Path to cortical neuron experimental Hi-C .npy")

    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    os.makedirs(args.output, exist_ok=True)

    if args.conditions:
        from configs.parameters import get_condition
        conditions = [get_condition(c) for c in args.conditions]
    else:
        conditions = SIMULATION_CONDITIONS

    if not args.analyze_only:
        logger.info(f"Running {len(conditions)} conditions × {args.n_replicates} replicates")
        logger.info(f"Engine: {args.engine}")
        run_all_simulations(conditions, args.n_replicates, args.engine, args.output, args.gpu)

    # Analyze
    logger.info("\nRunning analysis...")
    analyze_all_results(args.output, args.exp_es, args.exp_cn)

    # Plot
    logger.info("\nGenerating plots...")
    from analysis.plot_results import main as plot_main
    sys.argv = ["plot_results.py",
                "--results-dir", args.output,
                "--output", os.path.join(args.output, "figures")]
    if args.exp_es or args.exp_cn:
        data_dir = os.path.dirname(args.exp_es or args.exp_cn)
        sys.argv.extend(["--data-dir", data_dir])

    try:
        plot_main()
    except Exception as e:
        logger.warning(f"Plotting failed (matplotlib may not be available): {e}")


if __name__ == "__main__":
    main()
