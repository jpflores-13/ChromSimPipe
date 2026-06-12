#!/usr/bin/env python
"""
Merge shard outputs into a single trajectory per condition/replicate.

After multi-GPU parallel simulations, each shard produces its own
polychrom HDF5Reporter output (blocks_*.h5 files) or LEF-only contact maps.
This script merges them into a unified output directory.

For polychrom 3D conformations: copies all block files, renumbering frames
    so they form a continuous sequence across shards.
For LEF contact maps: sums the matrices (more samples = better statistics).

Usage:
    python merge_shards.py \
        --condition mESC_ctrl --replicate 3 --n-shards 4 \
        --results-dir results/polychrom_3d
"""

import os
import sys
import argparse
import logging
import json
import glob
import shutil
import numpy as np
import h5py
from multiprocessing import Pool, cpu_count

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import get_condition

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def merge_polychrom_hdf5(shard_dirs, output_dir):
    """
    Merge polychrom HDF5Reporter outputs (blocks_*.h5) from multiple shards.

    Each shard has its own blocks_0-99.h5, blocks_100-199.h5, etc.
    Instead of loading all frames into memory, we stream shard-by-shard
    and write immediately — keeping memory usage constant regardless of
    total frame count.
    """
    try:
        from polychrom.hdf5_format import list_URIs, load_URI, HDF5Reporter
        HAS_POLYCHROM = True
    except ImportError:
        HAS_POLYCHROM = False

    if not HAS_POLYCHROM:
        logger.warning("polychrom not installed — falling back to raw h5 merge")
        return _merge_raw_h5(shard_dirs, output_dir)

    # First pass: count total frames (cheap — only reads metadata)
    total_frames = 0
    shard_uri_counts = []
    for shard_dir in shard_dirs:
        try:
            uris = list_URIs(shard_dir)
            n = len(uris)
            shard_uri_counts.append((shard_dir, n))
            total_frames += n
            logger.info(f"  {os.path.basename(shard_dir)}: {n} frames")
        except Exception as e:
            logger.warning(f"  Failed to read {shard_dir}: {e}")
            shard_uri_counts.append((shard_dir, 0))

    if total_frames == 0:
        logger.warning("  No polychrom frames found to merge!")
        return

    # Stream: load one frame at a time, write immediately, then discard
    reporter = HDF5Reporter(folder=output_dir, max_data_length=100,
                            overwrite=True, blocks_only=True)

    block_idx = 0
    for shard_dir, n_uris in shard_uri_counts:
        if n_uris == 0:
            continue
        uris = list_URIs(shard_dir)
        for uri in uris:
            data = load_URI(uri)
            data["block"] = block_idx
            reporter.report("data", data)
            block_idx += 1

    reporter.dump_data()
    logger.info(f"  Merged {block_idx} total frames → {output_dir}")


def _merge_raw_h5(shard_dirs, output_dir):
    """
    Fallback: merge conformations.h5 files (old format) from multiple shards.

    Streams frames directly from source to destination HDF5 using h5py.copy()
    to avoid loading all frames into memory at once.
    """
    merged_path = os.path.join(output_dir, "conformations.h5")
    N = None
    total_frames = 0

    with h5py.File(merged_path, "w") as out_hf:
        for shard_dir in shard_dirs:
            h5_path = os.path.join(shard_dir, "conformations.h5")
            if not os.path.exists(h5_path):
                continue

            with h5py.File(h5_path, "r") as src_hf:
                n_frames = src_hf.attrs.get("n_frames", 0)
                if N is None:
                    N = src_hf.attrs.get("N", None)
                for i in range(n_frames):
                    src_key = f"frame_{i}"
                    if src_key in src_hf:
                        dst_key = f"frame_{total_frames}"
                        src_hf.copy(src_key, out_hf, name=dst_key)
                        total_frames += 1

            logger.info(f"  Copied {n_frames} frames from {os.path.basename(shard_dir)}")

        if total_frames == 0:
            logger.warning("  No conformation frames found to merge!")
            return

        out_hf.attrs["N"] = N or out_hf["frame_0"].shape[0]
        out_hf.attrs["n_frames"] = total_frames
        out_hf.attrs["merged_from_shards"] = len(shard_dirs)

    logger.info(f"  Merged {total_frames} total frames → {merged_path}")


def merge_lef_contact_maps(shard_dirs, output_dir):
    """
    Merge LEF-only contact maps by summing them.
    More shards = more samples = smoother map.
    """
    merged_map = None
    n_merged = 0

    for shard_dir in shard_dirs:
        npy_path = os.path.join(shard_dir, "lef_contact_map.npy")
        if not os.path.exists(npy_path):
            continue

        cmap = np.load(npy_path)
        if merged_map is None:
            merged_map = cmap.copy()
        else:
            merged_map += cmap
        n_merged += 1

    if merged_map is None:
        logger.warning("  No LEF contact maps found to merge!")
        return

    np.save(os.path.join(output_dir, "lef_contact_map.npy"), merged_map)
    logger.info(f"  Merged {n_merged} LEF contact maps")


def merge_params(shard_dirs, output_dir):
    """Copy params.json from the first shard, noting the merge."""
    for shard_dir in shard_dirs:
        params_path = os.path.join(shard_dir, "params.json")
        if os.path.exists(params_path):
            with open(params_path) as f:
                params = json.load(f)
            params["merged"] = True
            params["n_shards_merged"] = len(shard_dirs)
            with open(os.path.join(output_dir, "params.json"), "w") as f:
                json.dump(params, f, indent=2)
            return


def detect_output_format(shard_dirs):
    """Detect whether shards used polychrom HDF5Reporter or raw h5/npy."""
    for d in shard_dirs:
        if not os.path.exists(d):
            continue
        # Polychrom HDF5Reporter creates blocks_*.h5 files
        block_files = glob.glob(os.path.join(d, "blocks_*.h5"))
        if block_files:
            return "polychrom"
        # Old format: single conformations.h5
        if os.path.exists(os.path.join(d, "conformations.h5")):
            return "raw_h5"
        # LEF-only: contact map
        if os.path.exists(os.path.join(d, "lef_contact_map.npy")):
            return "lef_only"
    return "unknown"


def merge_one(results_dir, base_name, shard_dirs, cleanup=False):
    """Merge shards for a single condition/replicate.

    Parameters
    ----------
    results_dir : str
        Root results directory.
    base_name : str
        Condition/replicate base name (e.g. mESC_ctcf-mESC_rep0). The
        merged output directory is named ``merged_<base_name>`` so it is
        immediately distinguishable from individual shard directories.
    shard_dirs : list of str
        Full paths to all shard directories for this group.
    cleanup : bool
        If True, delete shard directories after a successful merge.
    """
    existing = [d for d in shard_dirs if os.path.exists(d)]
    if not existing:
        return

    output_dir = os.path.join(results_dir, f"merged_{base_name}")

    # Smart skip: only skip if we have already merged ALL current shards.
    # If new shard dirs have appeared since the last merge, re-merge everything.
    if os.path.exists(output_dir):
        merged_blocks = glob.glob(os.path.join(output_dir, "blocks_*.h5"))
        if merged_blocks:
            n_previously_merged = 0
            params_path = os.path.join(output_dir, "params.json")
            if os.path.exists(params_path):
                try:
                    with open(params_path) as f:
                        merged_params = json.load(f)
                    n_previously_merged = merged_params.get("n_shards_merged", 0)
                except Exception:
                    pass

            if len(existing) <= n_previously_merged:
                logger.info(f"  {base_name}: already merged "
                            f"({n_previously_merged} shards → {len(merged_blocks)} block files), "
                            f"skipping")
                return
            else:
                logger.info(f"  {base_name}: {len(existing)} shard dirs found but only "
                            f"{n_previously_merged} previously merged — re-merging all...")

    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"  {base_name}: merging {len(existing)} shards...")

    fmt = detect_output_format(existing)

    if fmt == "polychrom":
        merge_polychrom_hdf5(existing, output_dir)
    elif fmt == "raw_h5":
        _merge_raw_h5(existing, output_dir)
    elif fmt == "lef_only":
        merge_lef_contact_maps(existing, output_dir)
    else:
        merge_polychrom_hdf5(existing, output_dir)
        merge_lef_contact_maps(existing, output_dir)

    merge_params(existing, output_dir)

    if cleanup:
        for shard_dir in existing:
            shutil.rmtree(shard_dir)
            logger.info(f"  Deleted shard dir: {os.path.basename(shard_dir)}")


def _merge_one_wrapper(args):
    """Wrapper for multiprocessing Pool — unpacks args tuple."""
    results_dir, base_name, shard_dirs, cleanup = args
    try:
        merge_one(results_dir, base_name, shard_dirs, cleanup=cleanup)
    except Exception as e:
        logger.error(f"  Failed to merge {base_name}: {e}")


def discover_and_merge_all(results_dir, parallel=True, cleanup=False):
    """
    Auto-discover all shard directories in results_dir and merge them
    by condition/replicate.

    Shard dirs follow the naming pattern:
        <params>_ctcf-<type>_rep<N>_shard<M>

    Groups them by <params>_ctcf-<type>_rep<N> and merges each group.
    When parallel=True, merges multiple groups concurrently.
    """
    import re

    if not os.path.isdir(results_dir):
        logger.error(f"Results directory not found: {results_dir}")
        sys.exit(1)

    # Find all shard directories
    shard_pattern = re.compile(r'^(.+_rep\d+)_shard(\d+)$')
    groups = {}  # base_name → list of shard dirs

    for entry in sorted(os.listdir(results_dir)):
        full_path = os.path.join(results_dir, entry)
        if not os.path.isdir(full_path):
            continue
        m = shard_pattern.match(entry)
        if m:
            base_name = m.group(1)
            groups.setdefault(base_name, []).append(full_path)

    if not groups:
        logger.info("No shard directories found to merge.")
        return

    logger.info(f"Found {len(groups)} condition/replicate groups to merge:")
    for base_name, shard_dirs in sorted(groups.items()):
        logger.info(f"  {base_name}: {len(shard_dirs)} shards")

    # Sort shard dirs within each group
    work_items = []
    for base_name, shard_dirs in sorted(groups.items()):
        shard_dirs.sort()  # ensure shard0, shard1, shard2... order
        work_items.append((results_dir, base_name, shard_dirs, cleanup))

    if cleanup:
        logger.warning("--cleanup is set: shard directories will be deleted after merge!")

    if parallel and len(work_items) > 1:
        # Cap default = number of CPUs allocated to the job (SLURM_CPUS_PER_TASK
        # if set, falls back to multiprocessing.cpu_count()). Override with
        # MERGE_MAX_WORKERS=<n> in the environment.
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        cap_default = int(slurm_cpus) if slurm_cpus and slurm_cpus.isdigit() else cpu_count()
        cap = int(os.environ.get("MERGE_MAX_WORKERS", cap_default))
        n_workers = min(len(work_items), cap)
        logger.info(f"Merging {len(work_items)} groups in parallel ({n_workers} workers)")
        with Pool(n_workers) as pool:
            pool.map(_merge_one_wrapper, work_items)
    else:
        for item in work_items:
            _merge_one_wrapper(item)

    logger.info("All merges complete.")


def main():
    parser = argparse.ArgumentParser(
        description="Merge shard outputs",
        epilog="Examples:\n"
               "  # Merge ALL conditions/replicates found in results dir:\n"
               "  python merge_shards.py --all --results-dir results/polychrom_3d\n\n"
               "  # Merge a specific condition/replicate:\n"
               "  python merge_shards.py --condition mESC_ctrl --replicate 3 "
               "--n-shards 2 --results-dir results/polychrom_3d\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--all", action="store_true",
                        help="Auto-discover and merge ALL shard groups in results-dir")
    parser.add_argument("--condition", type=str, default=None)
    parser.add_argument("--replicate", type=int, default=None)
    parser.add_argument("--n-shards", type=int, default=None)
    parser.add_argument("--results-dir", type=str, required=True)
    parser.add_argument("--no-parallel", action="store_true",
                        help="Disable parallel merging (merge groups sequentially)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Delete shard directories after a successful merge. "
                             "Saves disk space but is irreversible — use with caution.")

    args = parser.parse_args()

    if args.all:
        discover_and_merge_all(args.results_dir,
                               parallel=not args.no_parallel,
                               cleanup=args.cleanup)
    elif args.condition and args.replicate is not None and args.n_shards:
        # Single condition/replicate mode (original behavior)
        cond = get_condition(args.condition)
        params_name = cond["params"]["name"]
        ctcf_type = cond["ctcf_type"]
        dir_prefix = f"{params_name}_ctcf-{ctcf_type}"

        shard_dirs = []
        for s in range(args.n_shards):
            d = os.path.join(args.results_dir,
                             f"{dir_prefix}_rep{args.replicate}_shard{s}")
            shard_dirs.append(d)

        existing = [d for d in shard_dirs if os.path.exists(d)]
        if not existing:
            logger.error(f"No shard directories found for "
                         f"{dir_prefix}_rep{args.replicate}")
            sys.exit(1)

        logger.info(f"Merging {len(existing)}/{args.n_shards} shards for "
                    f"{args.condition} ({dir_prefix}_rep{args.replicate})")

        merge_one(args.results_dir, f"{dir_prefix}_rep{args.replicate}", existing,
                  cleanup=args.cleanup)
    else:
        parser.error("Use --all, or provide --condition + --replicate + --n-shards")


if __name__ == "__main__":
    main()
