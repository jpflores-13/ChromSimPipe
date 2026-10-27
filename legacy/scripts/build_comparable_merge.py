#!/usr/bin/env python3
"""
build_comparable_merge.py — produce subsampled rep dirs for cross-condition
                            comparability (MSD CIs, fair pooled contact maps).

Why
---
After all conditions have been simulated, per-rep block counts vary because of
walltime and crash truncation, AND the number of reps per condition varies.
For (a) MSD per-rep alpha confidence intervals to be comparable across
conditions and (b) pooled contact maps to be averaged over equal samples,
we need every condition represented by the SAME number of replicates with
the SAME number of saved frames per replicate.

This script picks N reps per condition, takes a fixed F-block slice from
each rep's blocks_*.h5 files, and writes the slice as a new rep dir whose
condition name is the original suffixed with --suffix (default "_comparable").
The analysis pipeline then sees `<cond>_comparable` as its own condition
(parse_condition_rep just splits on the trailing `_repN`), and produces a
parallel set of `<cond>_comparable_<NNNN>blk_pooled_*` artifacts that you
can plot side-by-side with the full-pool ones for the original conditions.

Defaults (see CLAUDE.md "Statistical-depth policy"):
    --reps 3  --blocks-per-rep 700  --mode tail
That gives every condition 3 reps × 700 blocks × 100 frames = 210,000 pooled
frames (slight overshoot of the 200k cap; we snap to whole .h5 file
boundaries since polychrom packs 100 frames per file). Times 28 tiles
(TILING['n_tiles']) = 5.9M tile-conformations per condition ≈ 14× the
Gabriele 2022 Science-paper effective depth.

The 3-rep / 700-blocks-per-rep choice matches the catch-up policy
(MIN_REPS_PER_CONDITION=3 and BLOCKS_PER_REP=700 in
cluster/submit_catch_up.sh), so every catchup-completed rep contributes
exactly one tail-700-block slice to the comparable pool.

Storage
-------
By default the per-block .h5 files are SYMLINKED, not copied — each is
~150 MB and an 8-condition × 5-rep × 4-block subsample would otherwise be
~24 GB of duplicated data. Pass --copy if you need standalone files (e.g.
to ship the comparable subset to another machine without touching the
original sim dirs).

The polychrom HDF5Reporter writes 100 frames per blocks_*.h5 file; this
script's --blocks-per-rep is in those file units. Frame count = blocks * 100.

Selection mode
--------------
- tail (default): take the LAST n blocks of each chosen rep. Most robust
  for MSD comparability — reps may differ in their *first* block index due
  to warmup ordering, but the most-equilibrated frames sit at the tail.
- head: take the FIRST n blocks. Simpler indexing; use only if you know
  every rep's earliest blocks are already past the polymer-relaxation phase.

Usage
-----
    cd cohesin_sim
    python scripts/build_comparable_merge.py --dry-run
    python scripts/build_comparable_merge.py
    python scripts/build_comparable_merge.py --reps 3 --blocks-per-rep 1000

After it finishes, refresh the analysis layer:
    rm -f results/analysis/*_comparable_*
    sbatch cluster/submit_analysis_all_32cpu.sh
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

REP_DIR_RE = re.compile(r"^(?:merged_)?(.+)_rep(\d+)$")
SHARD_DIR_RE = re.compile(r"_rep\d+_shard\d+$")
BLOCK_INDEX_RE = re.compile(r"blocks_(\d+)")

# Mirrors DISK_TO_CONDITION in cluster/submit_catch_up.sh — these are the
# disk stems we recognise. Anything else (e.g. _common_smoke, legacy dirs)
# is ignored.
KNOWN_CONDITIONS = {
    "mESC_ctcf-mESC",
    "mESC_ctcf-neuron",
    "CN_baseline_ctcf-neuron",
    "CN_long_residency_ctcf-neuron",
    "CN_very_long_residency_ctcf-neuron",
    "CN_high_density_ctcf-neuron",
    "CN_long_res_high_dens_ctcf-neuron",
    "CN_long_res_low_dens_ctcf-neuron",
}


def discover_reps(results_dir: Path, suffix: str) -> dict[str, list[tuple[int, Path, int]]]:
    """Return {condition: [(rep_num, rep_dir, n_blocks), ...]} sorted by n_blocks desc.

    Skips shard dirs and any condition stem already ending in `suffix`
    (so re-running this tool doesn't recurse into its own outputs).

    When both `merged_<cond>_repN` and `<cond>_repN` exist for the same
    rep (legacy from before merge_shards.py:212 began prefixing with
    `merged_` on 2026-04-27), the merged_ form wins because it's the
    canonical post-refactor copy and bare-name leftovers are typically
    older partials with fewer blocks.
    """
    reps_by_cond: dict[str, list[tuple[int, Path, int]]] = defaultdict(list)
    seen: set[tuple[str, int]] = set()

    def _scan(d: Path) -> None:
        if not d.is_dir() or SHARD_DIR_RE.search(d.name):
            return
        m = REP_DIR_RE.match(d.name)
        if not m:
            return
        cond, rep = m.group(1), int(m.group(2))
        if cond.endswith(suffix) or cond not in KNOWN_CONDITIONS:
            return
        if (cond, rep) in seen:
            return
        seen.add((cond, rep))
        n_blocks = len(list(d.glob("blocks_*.h5")))
        if n_blocks > 0:
            reps_by_cond[cond].append((rep, d, n_blocks))

    # Pass 1: merged_* (canonical)
    for d in sorted(results_dir.glob("merged_*")):
        _scan(d)
    # Pass 2: bare-name (skipped if its merged_ counterpart was already seen)
    for d in sorted(results_dir.iterdir()):
        if d.name.startswith("merged_"):
            continue
        _scan(d)

    for cond in reps_by_cond:
        reps_by_cond[cond].sort(key=lambda t: (-t[2], t[0]))
    return reps_by_cond


def select_blocks(rep_dir: Path, n_blocks: int, mode: str) -> list[Path]:
    """Return up to n_blocks block files from rep_dir, ordered by block index."""
    files = sorted(
        rep_dir.glob("blocks_*.h5"),
        key=lambda p: int(BLOCK_INDEX_RE.search(p.name).group(1)),
    )
    if mode == "head":
        return files[:n_blocks]
    if mode == "tail":
        return files[-n_blocks:]
    raise ValueError(f"unknown mode: {mode!r}")


def materialise_rep(src_blocks: list[Path], src_dir: Path, dst_dir: Path,
                    use_symlinks: bool) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for src in src_blocks:
        dst = dst_dir / src.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if use_symlinks:
            os.symlink(os.path.relpath(src.resolve(), dst.parent), dst)
        else:
            shutil.copy2(src, dst)
    # Carry over params.json so analysis sees the same per-rep config.
    params_src = src_dir / "params.json"
    if params_src.exists():
        shutil.copy2(params_src, dst_dir / "params.json")


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--results-dir", type=Path, default=Path("results/polychrom_3d"),
                    help="root holding merged_<cond>_repN/ dirs (default: %(default)s)")
    ap.add_argument("--reps", type=int, default=3,
                    help="reps per condition in the comparable set (default: 3, matches MIN_REPS_PER_CONDITION)")
    ap.add_argument("--blocks-per-rep", type=int, default=700,
                    help="block files per rep; each holds 100 frames (default: 700 = 70 000 frames)")
    ap.add_argument("--mode", choices=["head", "tail"], default="tail",
                    help="which slice of each rep to take (default: tail = post-equilibration)")
    ap.add_argument("--copy", action="store_true",
                    help="copy block files instead of symlinking (uses ~150 MB per block)")
    ap.add_argument("--suffix", default="_comparable",
                    help="suffix appended to condition names (default: %(default)s)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print plan without writing")
    args = ap.parse_args()

    results_dir = args.results_dir
    if not results_dir.is_dir():
        print(f"ERROR: {results_dir} not found", file=sys.stderr)
        return 1

    reps_by_cond = discover_reps(results_dir, args.suffix)

    plan: list[dict] = []
    warnings: list[str] = []

    for cond in sorted(KNOWN_CONDITIONS):
        eligible = [
            t for t in reps_by_cond.get(cond, [])
            if t[2] >= args.blocks_per_rep
        ]
        chosen = eligible[: args.reps]

        if not chosen:
            avail = reps_by_cond.get(cond, [])
            avail_str = (
                ", ".join(f"rep{r}={n}blk" for r, _, n in avail) or "no reps"
            )
            warnings.append(
                f"  {cond}: no reps with >= {args.blocks_per_rep} blocks "
                f"(available: {avail_str}); SKIPPED"
            )
            continue
        if len(chosen) < args.reps:
            warnings.append(
                f"  {cond}: only {len(chosen)}/{args.reps} reps available with "
                f">= {args.blocks_per_rep} blocks — comparable set will be short."
            )

        for new_rep_idx, (orig_rep, src_dir, n_blocks) in enumerate(chosen, start=1):
            dst_name = f"merged_{cond}{args.suffix}_rep{new_rep_idx}"
            dst_dir = results_dir / dst_name
            blocks = select_blocks(src_dir, args.blocks_per_rep, args.mode)
            plan.append({
                "condition": cond,
                "comparable_rep": new_rep_idx,
                "src_rep_dir": str(src_dir),
                "src_orig_rep": orig_rep,
                "src_n_blocks": n_blocks,
                "blocks_taken": len(blocks),
                "dst_dir": str(dst_dir),
                "first_block": blocks[0].name if blocks else None,
                "last_block": blocks[-1].name if blocks else None,
            })

    print("=" * 72)
    print("build_comparable_merge plan")
    print("=" * 72)
    print(f"results_dir   : {results_dir}")
    print(f"reps/cond     : {args.reps}")
    print(f"blocks/rep    : {args.blocks_per_rep}  ({args.blocks_per_rep * 100} frames)")
    print(f"mode          : {args.mode}")
    print(f"link mode     : {'copy' if args.copy else 'symlink'}")
    print(f"output suffix : {args.suffix}")
    print(f"plan size     : {len(plan)} new rep dirs")
    print()

    if warnings:
        print("WARNINGS:")
        for w in warnings:
            print(w)
        print()

    by_cond: dict[str, list[dict]] = defaultdict(list)
    for p in plan:
        by_cond[p["condition"]].append(p)

    for cond in sorted(by_cond):
        ps = by_cond[cond]
        total_frames = len(ps) * args.blocks_per_rep * 100
        print(f"[{cond}{args.suffix}]  {len(ps)} reps -> {total_frames} pooled frames")
        for p in ps:
            print(
                f"  rep{p['comparable_rep']:<2} <- "
                f"{Path(p['src_rep_dir']).name} "
                f"(orig rep{p['src_orig_rep']}, {p['src_n_blocks']} blocks; "
                f"taking {p['blocks_taken']} {args.mode}: "
                f"{p['first_block']}..{p['last_block']})"
            )
        print()

    if args.dry_run:
        print("DRY RUN — not writing.")
        return 0
    if not plan:
        print("Nothing to do.")
        return 0

    for p in plan:
        src_dir = Path(p["src_rep_dir"])
        dst_dir = Path(p["dst_dir"])
        blocks = select_blocks(src_dir, args.blocks_per_rep, args.mode)
        materialise_rep(blocks, src_dir, dst_dir, use_symlinks=not args.copy)
        print(f"  materialised {dst_dir.name}  ({len(blocks)} blocks)")

    manifest_path = results_dir / "comparable_merge_manifest.json"
    manifest = {
        "results_dir": str(results_dir),
        "reps_per_condition": args.reps,
        "blocks_per_rep": args.blocks_per_rep,
        "frames_per_rep": args.blocks_per_rep * 100,
        "mode": args.mode,
        "link_mode": "copy" if args.copy else "symlink",
        "suffix": args.suffix,
        "warnings": warnings,
        "entries": plan,
    }
    with manifest_path.open("w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"\nmanifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
