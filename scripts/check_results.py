#!/usr/bin/env python3
"""
Check simulation results: list all completed replicates per condition.

Usage:
    python scripts/check_results.py
    python scripts/check_results.py --results-dir results/polychrom_3d
"""

import os
import sys
import argparse
import json
from collections import defaultdict

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from configs.parameters import SIMULATION_CONDITIONS


def check_results(results_dir):
    """Scan results directory and report what's available per condition."""

    if not os.path.isdir(results_dir):
        print(f"Results directory not found: {results_dir}")
        return

    # Collect all result directories
    all_dirs = sorted(os.listdir(results_dir))

    # Group by condition
    condition_names = [c["name"] for c in SIMULATION_CONDITIONS]
    results = defaultdict(list)   # condition_name -> list of (rep, dir_name, has_conformations)
    orphan_dirs = []

    for d in all_dirs:
        full_path = os.path.join(results_dir, d)
        if not os.path.isdir(full_path):
            continue

        # Skip shard directories (they should have been merged)
        if "_shard" in d:
            continue

        # Try to match to a condition
        matched = False
        for cond in SIMULATION_CONDITIONS:
            params_name = cond["params"]["name"]
            ctcf_type = cond["ctcf_type"]
            prefix = f"{params_name}_ctcf-{ctcf_type}_rep"

            if d.startswith(prefix):
                rep_str = d[len(prefix):]
                try:
                    rep = int(rep_str)
                except ValueError:
                    continue

                # Check if there are actual conformations
                has_conf = any(
                    f.endswith((".npy", ".h5", ".hdf5", ".xyz"))
                    for f in os.listdir(full_path)
                ) if os.listdir(full_path) else False

                # Check for params.json
                has_params = os.path.exists(os.path.join(full_path, "params.json"))

                results[cond["name"]].append({
                    "rep": rep,
                    "dir": d,
                    "has_conformations": has_conf,
                    "has_params": has_params,
                    "n_files": len(os.listdir(full_path)),
                })
                matched = True
                break

        if not matched:
            orphan_dirs.append(d)

    # Count leftover shard directories
    shard_dirs = [d for d in all_dirs if "_shard" in d and os.path.isdir(os.path.join(results_dir, d))]

    # --- Report ---
    print("=" * 70)
    print("SIMULATION RESULTS SUMMARY")
    print(f"Directory: {os.path.abspath(results_dir)}")
    print("=" * 70)
    print()

    total_complete = 0
    total_empty = 0

    for cond in SIMULATION_CONDITIONS:
        name = cond["name"]
        reps = sorted(results[name], key=lambda r: r["rep"])

        if reps:
            rep_nums = [r["rep"] for r in reps]
            complete = [r for r in reps if r["has_conformations"]]
            empty = [r for r in reps if not r["has_conformations"]]

            status = f"{len(complete)}/{len(reps)} complete"
            if empty:
                status += f" ({len(empty)} empty)"

            print(f"  {name}")
            print(f"    Replicates: {rep_nums}  —  {status}")
            for r in reps:
                marker = "OK" if r["has_conformations"] else "EMPTY"
                print(f"      rep{r['rep']:2d}: {r['n_files']:4d} files  [{marker}]")
            total_complete += len(complete)
            total_empty += len(empty)
        else:
            print(f"  {name}")
            print(f"    Replicates: NONE")
        print()

    print("-" * 70)
    print(f"Total: {total_complete} complete replicates, {total_empty} empty")

    if shard_dirs:
        print(f"\nLeftover shard directories ({len(shard_dirs)}):")
        for s in shard_dirs[:10]:
            print(f"  {s}")
        if len(shard_dirs) > 10:
            print(f"  ... and {len(shard_dirs) - 10} more")

    if orphan_dirs:
        print(f"\nUnrecognized directories ({len(orphan_dirs)}):")
        for o in orphan_dirs[:10]:
            print(f"  {o}")

    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check simulation results")
    parser.add_argument("--results-dir", default="results/polychrom_3d",
                        help="Results directory (default: results/polychrom_3d)")
    args = parser.parse_args()
    check_results(args.results_dir)
