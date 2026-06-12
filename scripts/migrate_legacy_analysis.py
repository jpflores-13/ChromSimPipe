#!/usr/bin/env python
"""Migrate per-sim-dir legacy analysis outputs to the flat results/analysis/ layout.

Before the 2026-04-27 flat-output refactor, run_analysis_all.py wrote each
replicate's heavy outputs to ``<sim_dir>/analysis/`` with un-prefixed names
(e.g. ``sim_contact_map.npy``, ``msd_overlay.png``) and copied prefixed
versions into ``results/analysis/``. After the refactor, every output goes
directly to ``results/analysis/{condition}_{n_blocks}blk_rep{rep}_*``.

This one-shot script lets you migrate older runs that still have data
under ``<sim_dir>/analysis/`` so it lives alongside the new outputs.

For each merged simulation directory under ``--results-dir``:
  1. compute the same prefix the live pipeline would (via ``_build_prefix``).
  2. for every file in ``<sim_dir>/analysis/``, move it to
     ``<analysis-dir>/{prefix}_{renamed}.{ext}`` using the same renaming
     rules ``_copy_to_common`` used (e.g. ``sim_contact_map.npy`` →
     ``{prefix}_contact_map.npy``, ``comparison_metrics.json`` →
     ``{prefix}_metrics.json``).
  3. optionally remove the now-empty ``<sim_dir>/analysis/`` directory.
  4. optionally rename ``<sim_dir>`` itself to ``merged_<sim_dir>`` so the
     new naming convention applies.

Safe to re-run: skips files whose destination already exists.

Usage examples
--------------
Dry-run (recommended first pass — logs what would happen, touches nothing)::

    python scripts/migrate_legacy_analysis.py \\
        --results-dir results/polychrom_3d --dry-run

Move the files in place and clean up empty directories::

    python scripts/migrate_legacy_analysis.py \\
        --results-dir results/polychrom_3d --delete-empty

Move + rename the merged sim dirs to the new ``merged_<basename>`` form::

    python scripts/migrate_legacy_analysis.py \\
        --results-dir results/polychrom_3d --delete-empty --rename-to-merged

Use ``--copy`` to leave the source files in place (useful for a careful
two-step audit before committing to the move).
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import sys

# Make scripts/ importable so we can reuse the live pipeline's helpers.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from run_analysis_all import (  # noqa: E402
    LEGACY_OUTPUT_NAMES,
    _build_prefix,
    is_merged_dir,
    parse_condition_rep,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


# Same renaming the old _copy_to_common helper applied. Anything not in this
# table is kept as ``{prefix}_{original_filename}`` (covers per-pair MSD
# artefacts whose suffix is the variable label, e.g. ``msd_<probe>.json``).
RENAME_MAP = {
    "sim_contact_map.npy":     "contact_map.npy",
    "sim_ps_curve.npz":        "ps_curve.npz",
    "sim_insulation.npy":      "insulation.npy",
    "comparison_metrics.json": "metrics.json",
}


def _dst_name(file_prefix: str, src_filename: str) -> str:
    """Compute the prefixed destination filename for one legacy file."""
    if src_filename in RENAME_MAP:
        return f"{file_prefix}_{RENAME_MAP[src_filename]}"
    # All other legacy files keep their original suffix; we just prepend the
    # prefix. Covers e.g. ``ps_metrics.json``, ``msd_overlay.png``,
    # ``rg_timecourse.npz``, ``contact_map_with_ctcf.png``,
    # ``apa_convergent.png``, ``calibration.json``, plus per-pair
    # ``msd_<label>.{json,npz}`` whose label is variable.
    return f"{file_prefix}_{src_filename}"


def migrate_one(
    sim_dir: str,
    analysis_dir: str,
    *,
    copy: bool = False,
    dry_run: bool = False,
    delete_empty: bool = False,
) -> int:
    """Migrate one merged sim dir's legacy ``analysis/`` files.

    Returns the number of files moved (or that would be moved, in dry-run).
    """
    legacy_dir = os.path.join(sim_dir, "analysis")
    if not os.path.isdir(legacy_dir):
        return 0

    file_prefix = _build_prefix(sim_dir)
    n_moved = 0

    if not dry_run:
        os.makedirs(analysis_dir, exist_ok=True)

    for fn in sorted(os.listdir(legacy_dir)):
        src = os.path.join(legacy_dir, fn)
        if not os.path.isfile(src):
            # Skip subdirectories (e.g. an old pooled_subdir if present).
            continue

        dst_name = _dst_name(file_prefix, fn)
        dst = os.path.join(analysis_dir, dst_name)

        if os.path.exists(dst):
            logger.warning(
                f"  [{file_prefix}] dst exists, skipping: {dst_name}"
            )
            continue

        if dry_run:
            logger.info(f"  [dry-run] {fn} -> {dst_name}")
            n_moved += 1
        elif copy:
            shutil.copy2(src, dst)
            logger.info(f"  copied : {fn} -> {dst_name}")
            n_moved += 1
        else:
            shutil.move(src, dst)
            logger.info(f"  moved  : {fn} -> {dst_name}")
            n_moved += 1

    # Optionally remove the now-empty legacy directory.
    if delete_empty and not dry_run and not copy:
        try:
            remaining = os.listdir(legacy_dir)
            if not remaining:
                os.rmdir(legacy_dir)
                logger.info(f"  removed empty {legacy_dir}")
            else:
                logger.info(
                    f"  NOT removing {legacy_dir} (still has "
                    f"{len(remaining)} non-file items)"
                )
        except OSError as e:
            logger.warning(f"  couldn't remove {legacy_dir}: {e}")

    return n_moved


def rename_to_merged(sim_dir: str, *, dry_run: bool = False) -> str | None:
    """Rename ``<sim_dir>`` -> ``merged_<sim_dir>`` if not already prefixed.

    Returns the new path on success, None if nothing to do (already
    prefixed or rename impossible).
    """
    parent = os.path.dirname(sim_dir)
    base = os.path.basename(sim_dir)
    if base.startswith("merged_"):
        return None
    new_path = os.path.join(parent, f"merged_{base}")
    if os.path.exists(new_path):
        logger.warning(
            f"  rename skipped: {new_path} already exists (manual cleanup needed)"
        )
        return None
    if dry_run:
        logger.info(f"  [dry-run] rename: {base} -> merged_{base}")
        return new_path
    os.rename(sim_dir, new_path)
    logger.info(f"  renamed: {base} -> merged_{base}")
    return new_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results-dir", required=True,
        help="Root holding merged simulation directories "
             "(e.g. results/polychrom_3d).",
    )
    parser.add_argument(
        "--analysis-dir", default=None,
        help="Destination folder. Default: <results-dir>/../analysis "
             "(matches run_analysis_all.py's default).",
    )
    parser.add_argument(
        "--copy", action="store_true",
        help="Copy files instead of moving them. Lets you keep originals "
             "until you're sure the migration looks correct.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Log what would happen without touching the filesystem.",
    )
    parser.add_argument(
        "--delete-empty", action="store_true",
        help="Remove each <sim_dir>/analysis/ once it has been emptied "
             "by a successful move. Ignored with --copy or --dry-run.",
    )
    parser.add_argument(
        "--rename-to-merged", action="store_true",
        help="After a successful migration of a sim dir, rename it from "
             "<sim_dir> to merged_<sim_dir> so the new naming convention "
             "applies. Skips dirs already prefixed.",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.results_dir):
        logger.error(f"--results-dir does not exist: {args.results_dir}")
        return 1

    if args.analysis_dir is None:
        args.analysis_dir = os.path.join(
            os.path.dirname(args.results_dir.rstrip("/\\")), "analysis"
        )
    logger.info(f"Source       : {args.results_dir}")
    logger.info(f"Destination  : {args.analysis_dir}")
    logger.info(f"Mode         : "
                f"{'dry-run' if args.dry_run else ('copy' if args.copy else 'move')}")

    n_dirs = 0
    n_files_total = 0
    for entry in sorted(os.listdir(args.results_dir)):
        sim_dir = os.path.join(args.results_dir, entry)
        if not os.path.isdir(sim_dir):
            continue
        if not is_merged_dir(sim_dir):
            # Skips shard dirs and any dir that doesn't look like a
            # merged simulation output (no blocks_*.h5 / conformations.h5
            # / lef_contact_map.npy).
            continue
        if not os.path.isdir(os.path.join(sim_dir, "analysis")):
            continue

        condition, rep = parse_condition_rep(
            entry[len("merged_"):] if entry.startswith("merged_") else entry
        )
        logger.info(f"--- {entry}  (condition={condition}, rep={rep}) ---")

        moved = migrate_one(
            sim_dir, args.analysis_dir,
            copy=args.copy, dry_run=args.dry_run,
            delete_empty=args.delete_empty,
        )
        n_dirs += 1
        n_files_total += moved

        if args.rename_to_merged and not args.copy:
            rename_to_merged(sim_dir, dry_run=args.dry_run)

    logger.info("")
    logger.info(
        f"Done: processed {n_dirs} merged dir(s), "
        f"{'would move' if args.dry_run else ('copied' if args.copy else 'moved')} "
        f"{n_files_total} file(s) in total."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
