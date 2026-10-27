#!/usr/bin/env python
"""
Modular analysis pipeline for cohesin-simulation results.

Discovers every merged result directory (blocks_*.h5 / conformations.h5 /
lef_contact_map.npy, excluding shard dirs), then runs the full suite of
analyses per replicate and pooled per condition:

  - contact map                   (analysis/contact_maps.py)
  - P(s) curve + power-law fit    (analysis/ps_curve.py)
  - insulation score              (analysis/contact_maps.py)
  - absolute loop quantification  (analysis/absolute_quant.py)
  - CTCF relative BED + figure    (analysis/ctcf_plotting.py)
  - experimental comparison       (analysis/experimental_compare.py)
    Uses .npy legacy matrices *or* .mcool files, whichever is found.

Per-replicate outputs (in <sim_dir>/analysis/ and the common output dir):

    sim_contact_map.npy
    sim_ps_curve.npz
    sim_insulation.npy
    ps_metrics.json
    ctcf_sites_relative.bed
    contact_map_with_ctcf.png
    apa_convergent.png
    apa_loops_quant.json
    comparison_metrics.json          (if experimental matrix available)
    sim_vs_exp_map.png               (if mcool match found)

Pooled-per-condition outputs:

    <cond>_<blocks>blk_pooled_contact_map.npy
    <cond>_<blocks>blk_pooled_ps_curve.npz
    <cond>_<blocks>blk_pooled_insulation.npy
    <cond>_<blocks>blk_pooled_ctcf_sites_relative.bed
    <cond>_<blocks>blk_pooled_contact_map_with_ctcf.png
    <cond>_<blocks>blk_pooled_apa_convergent.png
    <cond>_<blocks>blk_pooled_ps_metrics.json
    <cond>_<blocks>blk_pooled_metrics.json        (if exp available)

Overlay plot across all conditions:

    ps_overlay_all_conditions.png

Usage
-----
# Analyse everything found under results/polychrom_3d/:
    python scripts/run_analysis_all.py --results-dir results/polychrom_3d

# Each condition is compared with its compare_hic map from the config,
# looked up as <hic-dir>/<hic_type>.cool (or .npy), else data/mcool/:
    python scripts/run_analysis_all.py \\
        --results-dir results/polychrom_3d \\
        --hic-dir data/mcool

Options
-------
--output-dir DIR     Common output folder (default: <results_dir>/../analysis)
--n-jobs N           Parallel workers for contact detection inside each map
--no-parallel        Analyse directories sequentially (debugging)
--skip-existing      Skip dirs that already have analysis/sim_contact_map.npy
--no-pool            Skip the pooled-replicate step
--no-apa             Skip absolute-loop-quant APA pileups
--no-ctcf-overlay    Skip the CTCF-aligned figure and relative BED
--mcool HIC_TYPE=PATH      Override the mcool used for a compare_hic type (repeatable)
--ctcf-bed CTCF_TYPE=PATH  Override the oriented CTCF BED for a ctcf_type (repeatable)
--elements-bed CTCF_TYPE=PATH
                           Optional BED of non-CTCF sticky elements (enhancers /
                           promoters) to overlay for that ctcf_type (repeatable).
--elements-label TEXT      Legend label for the sticky-element overlay.
                           Default: 'enhancers/promoters'.
"""

import os
import sys
import re
import glob
import json
import shutil
import argparse
import logging
import traceback
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import cpu_count


def _allocated_cpu_count() -> int:
    """CPUs actually available to this process.

    Order of trust:
      1. ``SLURM_CPUS_PER_TASK`` when running inside SLURM. This is the only
         value the scheduler *intends* us to use. On clusters whose
         TaskPlugin is not ``task/cgroup``, ``sched_getaffinity`` can
         include SMT siblings (32 logical for 16 physical), which silently
         doubles the Phase 1 ThreadPoolExecutor width — observed on
         2026-04-23 job 9388487 under --cpus-per-task=16 +
         --hint=nomultithread, where the orchestrator spawned 2
         concurrent dirs and re-introduced the parent-heap-fragmentation
         OOM path the sequential layout was meant to avoid.
      2. Linux cgroup/cpuset via ``sched_getaffinity`` for non-SLURM runs.
      3. ``multiprocessing.cpu_count()`` as a last resort.
    """
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus and slurm_cpus.isdigit():
        return int(slurm_cpus)
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return cpu_count()


# Allow imports from project root and analysis/
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "analysis"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Condition → CTCF type / Hi-C reference lookup (from config/ChromSimConfig.yaml)
# ---------------------------------------------------------------------------
# Merged directories are named merged_{params_name}_ctcf-{ctcf_type}_rep{N}
# (see scripts/merge_shards.py). The matching condition in the config says
# which CTCF set was simulated and which Hi-C map (compare_hic) to compare to.

from configs.parameters import (  # noqa: E402
    SIMULATION_CONDITIONS, CTCF_TYPES, CTCF_BEDS, MCOOL_FILES,
)

_CONDITION_BY_KEY = {
    f"{c['params']['name']}_ctcf-{c['ctcf_type']}": c for c in SIMULATION_CONDITIONS
}
_CTCF_TYPE_RE = re.compile(r"_ctcf-([^_]+)")


def parse_type_overrides(values: list[str] | None, flag: str) -> dict[str, str]:
    """Parse repeated TYPE=PATH command-line values into a dict."""
    out: dict[str, str] = {}
    for v in values or []:
        if "=" not in v:
            raise SystemExit(f"{flag} expects TYPE=PATH, got {v!r}")
        k, path = v.split("=", 1)
        out[k] = path
    return out


def lookup_condition(base_name: str) -> dict | None:
    """Return the config condition for a merged/pooled dir or condition name."""
    b = base_name[len("merged_"):] if base_name.startswith("merged_") else base_name
    b = re.sub(r"_rep\d+$", "", b)
    if b in _CONDITION_BY_KEY:
        return _CONDITION_BY_KEY[b]
    for key in sorted(_CONDITION_BY_KEY, key=len, reverse=True):
        if b.startswith(key):
            return _CONDITION_BY_KEY[key]
    return None


def infer_cell_type(base_name: str) -> str | None:
    """Return the ctcf_type a directory was simulated with, or None."""
    cond = lookup_condition(base_name)
    if cond:
        return cond["ctcf_type"]
    m = _CTCF_TYPE_RE.search(base_name)
    return m.group(1) if m else None


def find_hic_path(base_name: str, hic_dir: str | None) -> str | None:
    """Return an experimental Hi-C .npy (<hic_dir>/<compare_hic>.npy), or None."""
    cond = lookup_condition(base_name)
    if not hic_dir or not cond or not cond.get("compare_hic"):
        return None
    candidate = os.path.join(hic_dir, f"{cond['compare_hic']}.npy")
    return candidate if os.path.exists(candidate) else None


def default_ctcf_bed(cell_type: str | None, repo_root: str) -> str | None:
    """Return the oriented CTCF BED for a ctcf_type if it exists on disk."""
    rel = CTCF_BEDS.get(cell_type) if cell_type else None
    if not rel:
        return None
    path = rel if os.path.isabs(rel) else os.path.join(repo_root, rel)
    return path if os.path.exists(path) else None


def resolve_ctcf_bed(
    base_name: str,
    overrides: dict[str, str],
    repo_root: str,
) -> tuple[str | None, str | None]:
    """(ctcf_type, bed_path); bed_path is None if no file was found."""
    ct = infer_cell_type(base_name)
    if ct in overrides:
        return ct, overrides[ct]
    return ct, default_ctcf_bed(ct, repo_root)


def find_mcool(hic_type: str, hic_dir: str | None,
               overrides: dict[str, str]) -> str | None:
    """Return the first existing mcool for a Hi-C type: override, hic_dir, data/mcool."""
    candidates = [overrides.get(hic_type)]
    if hic_dir:
        candidates.append(os.path.join(hic_dir, f"{hic_type}.cool"))
    candidates.append(os.path.join(_REPO_ROOT, MCOOL_FILES[hic_type]))
    return next((c for c in candidates if c and os.path.exists(c)), None)


def resolve_mcool(
    base_name: str,
    hic_dir: str | None,
    overrides: dict[str, str],
) -> str | None:
    """Return the mcool for the condition's compare_hic map, or None."""
    cond = lookup_condition(base_name)
    hic = cond.get("compare_hic") if cond else None
    if not hic:
        logger.info(f"    [HiC routing] {base_name!r}: no compare_hic in config")
        return None
    path = find_mcool(hic, hic_dir, overrides)
    if path:
        logger.info(f"    [HiC routing] {base_name!r} → {hic} ({path})")
    else:
        logger.warning(f"    [HiC routing] {base_name!r}: mcool for '{hic}' not found")
    return path


def resolve_elements_bed(base_name: str, overrides: dict[str, str]) -> str | None:
    """
    Return the "sticky-element" BED (enhancers / promoters / other) for the
    condition's ctcf_type, or None if the user didn't supply one.
    """
    return overrides.get(infer_cell_type(base_name))


# ---------------------------------------------------------------------------
# Directory name parsing
# ---------------------------------------------------------------------------

REP_RE = re.compile(r"^(.+)_rep(\d+)$")
MERGED_PREFIX = "merged_"


def parse_condition_rep(base_name: str) -> tuple[str, int]:
    """Parse 'ctcf_only_ctcf-condition2_rep3' → ('ctcf_only_ctcf-condition2', 3).

    A leading 'merged_' (the prefix written by scripts/merge_shards.py for
    merged simulation directories) is stripped before parsing so condition
    names — and therefore every analysis output filename — stay clean.
    """
    if base_name.startswith(MERGED_PREFIX):
        base_name = base_name[len(MERGED_PREFIX):]
    m = REP_RE.match(base_name)
    if not m:
        return base_name, -1
    return m.group(1), int(m.group(2))


def count_conformations(sim_dir: str) -> int:
    """Count the number of conformations in a merged dir without loading them all."""
    block_files = sorted(glob.glob(os.path.join(sim_dir, "blocks_*.h5")))
    if block_files:
        try:
            from polychrom.hdf5_format import list_URIs
            return len(list_URIs(sim_dir))
        except ImportError:
            pass
        # Manual count: read each block file's keys
        import h5py
        count = 0
        for bf in block_files:
            with h5py.File(bf, "r") as hf:
                count += len([k for k in hf.keys() if isinstance(hf[k], h5py.Group)])
        return count

    h5_path = os.path.join(sim_dir, "conformations.h5")
    if os.path.exists(h5_path):
        import h5py
        with h5py.File(h5_path, "r") as hf:
            return int(hf.attrs.get("n_frames", len(hf.keys())))

    return 0


def _build_prefix(sim_dir: str) -> str:
    """Per-replicate filename prefix for a merged simulation directory.

    Returns ``{condition}_{n_blocks}blk_rep{rep}`` when the directory
    name contains a ``_repN`` suffix, otherwise ``{base}_{n_blocks}blk``.
    Used by every analysis output so files in ``results/analysis/``
    sort by condition.
    """
    base = os.path.basename(sim_dir.rstrip("/\\"))
    condition, rep = parse_condition_rep(base)
    n_blocks = count_conformations(sim_dir)
    if rep < 0:
        return f"{base}_{n_blocks}blk"
    return f"{condition}_{n_blocks}blk_rep{rep}"


LEGACY_OUTPUT_NAMES = frozenset({
    "sim_contact_map.npy",
    "sim_ps_curve.npz",
    "sim_insulation.npy",
    "ps_metrics.json",
    "comparison_metrics.json",
    "ctcf_sites_relative.bed",
    "contact_map_with_ctcf.png",
    "apa_convergent.png",
    "apa_loops_quant.json",
    "sim_vs_exp_map.png",
    "msd_overlay.png",
    "rg_timecourse.png",
    "rg_timecourse.npz",
    "dwell_times.png",
    "pair_distance.png",
    "loop_fractions.json",
    "calibration.json",
})


def _warn_legacy_outputs(common_dir: str) -> None:
    """Warn if any un-prefixed legacy filenames are present in common_dir.

    These would indicate a regression where an analysis module wrote a file
    without applying the {condition}_*blk_*_ prefix.
    """
    if not os.path.isdir(common_dir):
        return
    for fn in os.listdir(common_dir):
        if fn in LEGACY_OUTPUT_NAMES:
            logger.warning(f"  [layout] un-prefixed legacy file in {common_dir}: {fn} "
                           f"(should have a condition prefix)")


# ---------------------------------------------------------------------------
# Directory discovery
# ---------------------------------------------------------------------------

SHARD_RE = re.compile(r"_shard\d+$")


# ---------------------------------------------------------------------------
# Two-point MSD + polymer dynamics runner
# ---------------------------------------------------------------------------

def _run_msd_and_dynamics(
    conformations,
    cell_type: str,
    common_dir: str,
    file_prefix: str,
    display_name: str = "",
    *,
    do_msd: bool = True,
    do_polymer_dynamics: bool = True,
    calibrate_with: str = "hic",
    expt_msd_csv: str | None = None,
    expt_ps_curve: tuple | None = None,
    sim_ps_curve: tuple | None = None,
) -> dict:
    """
    Run two-point MSD, polymer-dynamics diagnostics, and physical-units
    calibration on one conformation list (per-replicate OR pooled).

    Outputs (all in ``common_dir``, prefixed with ``file_prefix``):

        {prefix}_msd_<label>.json     MSD summary per pair
        {prefix}_msd_<label>.npz      MSD arrays per pair
        {prefix}_msd_overlay.png      overlay of all MSD pairs
        {prefix}_rg_timecourse.png    R_g(t) for the locus
        {prefix}_rg_timecourse.npz    raw R_g(t)
        {prefix}_loop_fractions.json  looped-fraction scalar per pair
        {prefix}_dwell_times.png      dwell-time histogram across pairs
        {prefix}_pair_distance.png    pair-separation density across pairs
        {prefix}_calibration.json     nm/monomer + sec/frame + source

    ``display_name`` is a short human label used ONLY in plot titles /
    log lines so it is clear which condition this block refers to.

    Parameters
    ----------
    conformations : list of (N_beads, 3) arrays OR None
        Short-circuits (returns empty dict) if None; LEF-only data cannot
        produce MSD curves.
    cell_type : str
        Used to pick MSD pairs and the matching experimental reference.
    common_dir : str
        Flat output folder (results/analysis/) where every artefact lands.
    file_prefix : str
        Stem prepended to every filename (e.g. ``cond1_baseline_2400blk_rep0`` or
        ``cond1_baseline_7200blk_pooled``).
    display_name : str, optional
        Short human label used only in plot titles and log lines.
    calibrate_with : {"hic", "msd", "none"}
        Physical-units anchor. "hic" uses sim/expt P(s); "msd" needs
        ``expt_msd_csv``; "none" skips calibration.
    expt_ps_curve, sim_ps_curve : (distances, ps) tuples, optional
        Needed for the "hic" calibration anchor. Pass None to fall back to
        an assumed nm-per-monomer from the polymer physics heuristic.

    Returns
    -------
    dict with the MSD curves per pair, the calibration dataclass, and
    per-pair metadata. Mostly useful for the pooled overlay across
    conditions.
    """
    import numpy as np
    if conformations is None:
        _label = display_name or cell_type
        logger.info(f"    [MSD skipped: no 3D conformations for {_label}]")
        return {"msd_curves": {}, "calibration": None,
                "pairs": [], "cell_type": cell_type}

    import msd_two_point as msd_mod
    import polymer_dynamics as pd_mod
    import calibration as cal_mod
    from configs.parameters import (
        RESOLUTION, TILING, N_MONOMERS, SIM_RUN,
        get_msd_pairs, MSD_FIT, CALIBRATION,
    )

    tile_size = TILING["tile_size"]
    pad = TILING["padding"]
    n_tiles = TILING["n_tiles"]
    contact_radius = SIM_RUN["contact_radius"]

    os.makedirs(common_dir, exist_ok=True)
    out: dict = {
        "msd_curves": {},
        "pairs": [],
        "cell_type": cell_type,
        "calibration": None,
    }

    pairs, labels = get_msd_pairs(cell_type)
    logger.info(f"    [MSD] {cell_type}: {len(pairs)} pair(s) "
                f"{list(zip(labels, pairs))}")

    # ── Two-point MSD per pair ──────────────────────────────────────────
    if do_msd:
        for (idx_a, idx_b), label in zip(pairs, labels):
            try:
                res = msd_mod.run_msd_for_pair(
                    conformations, pair=(idx_a, idx_b), label=label,
                    out_dir=common_dir, file_prefix=file_prefix,
                    tile_size=tile_size, pad=pad, n_tiles=n_tiles,
                    lag_min=MSD_FIT["lag_min_frames"],
                    lag_max_frac=MSD_FIT["lag_max_frac"],
                    fit_lag_min=MSD_FIT["fit_lag_min"],
                    fit_lag_max_frac=MSD_FIT["fit_lag_max_frac"],
                    min_n_lags_for_fit=MSD_FIT["min_n_lags_for_fit"],
                )
                out["msd_curves"][label] = res
                out["pairs"].append({"label": label, "a": int(idx_a),
                                     "b": int(idx_b),
                                     "sep_monomers": int(abs(idx_b - idx_a))})
            except Exception as e:  # noqa: BLE001
                logger.warning(f"    [MSD] pair {label} failed: {e}")

        if out["msd_curves"]:
            try:
                msd_mod.plot_msd_curves(
                    out["msd_curves"],
                    os.path.join(common_dir, f"{file_prefix}_msd_overlay.png"),
                    title=f"Two-point MSD: {display_name or cell_type}",
                )
            except Exception as e:  # noqa: BLE001
                logger.warning(f"    [MSD] overlay failed: {e}")

    # ── Polymer-dynamics diagnostics per pair ────────────────────────────
    if do_polymer_dynamics:
        try:
            rg = pd_mod.compute_rg_timecourse(
                conformations, bead_range=(0, N_MONOMERS),
                tile_size=tile_size, pad=pad, n_tiles=n_tiles,
            )
            pd_mod.plot_rg_timecourse(
                {display_name or cell_type: rg},
                os.path.join(common_dir, f"{file_prefix}_rg_timecourse.png"),
                title=f"Radius of gyration: {display_name or cell_type}",
            )
            np.savez_compressed(
                os.path.join(common_dir, f"{file_prefix}_rg_timecourse.npz"),
                t=rg["t"], rg=rg["rg"],
            )
            out["rg_mean"] = float(rg["rg"].mean())
            out["rg_std"] = float(rg["rg"].std())
        except Exception as e:  # noqa: BLE001
            logger.warning(f"    [R_g] failed: {e}")

        dwell_by_label: dict = {}
        pair_dist_by_label: dict = {}
        loop_fracs: list[dict] = []
        for (idx_a, idx_b), label in zip(pairs, labels):
            try:
                lf = pd_mod.compute_looped_fraction(
                    conformations, idx_a, idx_b,
                    contact_radius=contact_radius,
                    tile_size=tile_size, pad=pad, n_tiles=n_tiles,
                )
                loop_fracs.append({"label": label, **lf})

                dw = pd_mod.compute_dwell_times(
                    conformations, idx_a, idx_b,
                    contact_radius=contact_radius,
                    tile_size=tile_size, pad=pad, n_tiles=n_tiles,
                )
                dwell_by_label[label] = dw

                pdist = pd_mod.compute_pair_distance_distribution(
                    conformations, idx_a, idx_b,
                    tile_size=tile_size, pad=pad, n_tiles=n_tiles,
                )
                pair_dist_by_label[label] = pdist
            except Exception as e:  # noqa: BLE001
                logger.warning(f"    [dyn] pair {label} failed: {e}")

        out["loop_fracs"] = loop_fracs
        if loop_fracs:
            with open(os.path.join(common_dir, f"{file_prefix}_loop_fractions.json"), "w") as f:
                json.dump(loop_fracs, f, indent=2)
        if dwell_by_label:
            try:
                pd_mod.plot_dwell_time_hist(
                    dwell_by_label,
                    os.path.join(common_dir, f"{file_prefix}_dwell_times.png"),
                    title=f"Contact dwell times: {display_name or cell_type}",
                )
            except Exception as e:  # noqa: BLE001
                logger.warning(f"    [dwell] plot failed: {e}")
        if pair_dist_by_label:
            try:
                pd_mod.plot_pair_distance_hist(
                    pair_dist_by_label,
                    os.path.join(common_dir, f"{file_prefix}_pair_distance.png"),
                    title=f"Pair separation density: {display_name or cell_type}",
                )
            except Exception as e:  # noqa: BLE001
                logger.warning(f"    [pair-dist] plot failed: {e}")

    # ── Physical-units calibration ──────────────────────────────────────
    calib = None
    try:
        hic_cfg = CALIBRATION["hic"]
        if calibrate_with == "hic":
            if sim_ps_curve is not None and expt_ps_curve is not None:
                sim_d, sim_p = sim_ps_curve
                exp_d, exp_p = expt_ps_curve
                calib = cal_mod.calibrate_from_hic(
                    sim_d, sim_p, exp_d, exp_p,
                    s_ref_bp=hic_cfg["s_ref_bp"],
                    bp_per_monomer=hic_cfg["fixed_bp_per_monomer"],
                    persistence_nm=hic_cfg["persistence_nm"],
                )
            else:
                calib = cal_mod.get_calibration(
                    "assumed",
                    nm_per_monomer=hic_cfg["persistence_nm"]
                        * (hic_cfg["fixed_bp_per_monomer"] / 1000.0) ** 0.5,
                    sec_per_frame=1.0,
                )
                calib.source = "assumed (no expt P(s))"
        elif calibrate_with == "msd":
            if not expt_msd_csv or not os.path.exists(expt_msd_csv):
                logger.info(f"    [calibration] no expt MSD CSV ({expt_msd_csv!r}); "
                            f"skipping msd anchor")
            elif not out["msd_curves"]:
                logger.info("    [calibration] no sim MSD curve; skipping msd anchor")
            else:
                sim_lags, sim_msd = next(iter(out["msd_curves"].values()))["lags"], \
                                    next(iter(out["msd_curves"].values()))["msd"]
                exp_t, exp_m = cal_mod.load_experimental_msd_csv(expt_msd_csv)
                calib = cal_mod.calibrate_from_msd(
                    sim_lags, sim_msd, exp_t, exp_m,
                )
        elif calibrate_with == "none":
            pass
        else:
            logger.warning(f"    [calibration] unknown anchor {calibrate_with!r}; "
                           f"skipping")
    except Exception as e:  # noqa: BLE001
        logger.warning(f"    [calibration] failed: {e}")
        calib = None

    if calib is not None:
        out["calibration"] = calib
        with open(os.path.join(common_dir, f"{file_prefix}_calibration.json"), "w") as f:
            json.dump(calib.to_dict(), f, indent=2)
        logger.info(f"    [calibration] nm/monomer={calib.nm_per_monomer:.2f} "
                    f"sec/frame={calib.sec_per_frame:.4f} "
                    f"source={calib.source}")

    return out


def is_merged_dir(path: str) -> bool:
    """True if this directory looks like a merged (not shard) result dir."""
    name = os.path.basename(path)
    if SHARD_RE.search(name):
        return False            # it's a shard dir
    if not os.path.isdir(path):
        return False
    has_blocks = bool(glob.glob(os.path.join(path, "blocks_*.h5")))
    has_conf   = os.path.exists(os.path.join(path, "conformations.h5"))
    has_lef    = os.path.exists(os.path.join(path, "lef_contact_map.npy"))
    return has_blocks or has_conf or has_lef


def discover_merged_dirs(results_dir: str) -> list[str]:
    """Return sorted list of merged result directories inside results_dir."""
    merged = []
    for entry in sorted(os.listdir(results_dir)):
        full = os.path.join(results_dir, entry)
        if is_merged_dir(full):
            merged.append(full)
    return merged


# ---------------------------------------------------------------------------
# Per-replicate analysis runner
# ---------------------------------------------------------------------------

def run_analysis_one(sim_dir: str, hic_path: str | None,
                     n_jobs: int, skip_existing: bool,
                     common_dir: str | None = None,
                     ctcf_bed_path: str | None = None,
                     cell_type: str | None = None,
                     mcool_path: str | None = None,
                     do_apa: bool = True,
                     do_ctcf_overlay: bool = True,
                     elements_bed_path: str | None = None,
                     elements_label: str = "enhancers/promoters",
                     do_msd: bool = True,
                     do_polymer_dynamics: bool = True,
                     calibrate_with: str = "hic",
                     expt_msd_csv: str | None = None,
                     reuse_heavy: bool = False,
                     cleanup_legacy_cache: bool = False) -> bool:
    """
    Run the full modular analysis pipeline for one merged directory.

    Orchestrates: contact_maps (map + insulation), ps_curve (P(s) fit),
    absolute_quant (APA + per-loop AbLE), ctcf_plotting (relative BED +
    aligned figure), experimental_compare (legacy .npy OR .mcool).

    ``reuse_heavy`` enables "resume mode": if the CPU-heavy cached outputs
    (``sim_contact_map.npy``, ``rg_timecourse.npz`` + ``msd_overlay.png``)
    already exist in ``<sim_dir>/analysis``, they are loaded from disk
    instead of being recomputed from conformations. The cheap, downstream
    CTCF-dependent steps (APA pileup, CTCF overlay, experimental comparison)
    are *always* re-run so they pick up any new CTCF BED or settings.

    Returns True on success, False on failure.
    """
    import numpy as np
    from contact_maps import (
        load_conformations_h5, load_lef_contact_map,
        extract_tiles_and_average,
        compute_insulation_score,
        compare_contact_maps,
    )
    import ps_curve as ps_mod
    import absolute_quant as abs_mod
    import ctcf_plotting as ctcf_mod
    import experimental_compare as exp_mod

    from configs.parameters import RESOLUTION, CHROM, REGION_START, REGION_END

    base_name = os.path.basename(sim_dir)

    if common_dir is None:
        logger.error(f"  FAILED: {base_name}: common_dir is required (flat output layout).")
        return False

    file_prefix = _build_prefix(sim_dir)
    contact_map_path = os.path.join(common_dir, f"{file_prefix}_contact_map.npy")

    legacy_dir = os.path.join(sim_dir, "analysis")
    if os.path.isdir(legacy_dir):
        logger.warning(f"  [legacy cache] {legacy_dir} exists from a previous "
                       f"layout; it is no longer read or written. "
                       f"Re-run with CLEANUP_LEGACY=1 to remove it.")

    if skip_existing and os.path.exists(contact_map_path):
        logger.info(f"  [skip] {base_name}: {file_prefix}_contact_map.npy already exists")
        return True

    logger.info(f"  Analysing: {base_name}  (cell_type={cell_type})")
    if reuse_heavy:
        logger.info(f"    [resume mode] cached heavy outputs will be reused "
                    f"if present; downstream steps always re-run.")
    if hic_path:
        logger.info(f"    legacy Hi-C npy: {os.path.basename(hic_path)}")
    if mcool_path:
        logger.info(f"    mcool: {os.path.basename(mcool_path)}")
    if ctcf_bed_path:
        logger.info(f"    CTCF BED: {os.path.basename(ctcf_bed_path)}")
    if elements_bed_path:
        logger.info(f"    elements BED: {os.path.basename(elements_bed_path)} "
                    f"({elements_label})")

    try:
        os.makedirs(common_dir, exist_ok=True)

        # --- Load conformations → contact map ---
        # In resume mode, skip the expensive cKDTree contact extraction when
        # sim_contact_map.npy already exists on disk. Conformations are loaded
        # lazily later only if MSD / polymer-dynamics also need to be redone.
        conformations = None
        sim_map = None

        if reuse_heavy and os.path.exists(contact_map_path):
            sim_map = np.load(contact_map_path)
            logger.info(f"    {base_name}: [resume] loaded cached contact map "
                        f"({sim_map.shape[0]}×{sim_map.shape[1]})")
        else:
            try:
                conformations = load_conformations_h5(sim_dir)
            except FileNotFoundError:
                pass

            if conformations is not None:
                from configs.parameters import SIM_RUN
                contact_radius = SIM_RUN["contact_radius"]
                n = len(conformations) if hasattr(conformations, "__len__") else "?"
                logger.info(f"    {base_name}: {n} conformations, {n_jobs} workers")
                sim_map = extract_tiles_and_average(conformations, contact_radius, n_jobs)
            else:
                sim_map = load_lef_contact_map(sim_dir)
                if sim_map is None:
                    logger.error(f"  FAILED: {base_name}: no simulation data found")
                    return False
                logger.info(f"    {base_name}: using LEF bridging contact map")

            np.save(contact_map_path, sim_map)

        # --- P(s) via ps_curve module ---
        distances, ps = ps_mod.compute_ps_curve(sim_map, RESOLUTION)
        np.savez(os.path.join(common_dir, f"{file_prefix}_ps_curve.npz"),
                 distances=distances, ps=ps)
        ps_payload = ps_mod.save_ps_json(
            os.path.join(common_dir, f"{file_prefix}_ps_metrics.json"),
            distances, ps, RESOLUTION,
        )
        exp_alpha = ps_payload["metrics"].get("ps_exponent_2_30kb")
        if exp_alpha is not None:
            logger.info(f"    {base_name}: P(s) exponent α = {exp_alpha:.3f}")

        # --- Insulation score ---
        insulation = compute_insulation_score(sim_map)
        np.save(os.path.join(common_dir, f"{file_prefix}_insulation.npy"), insulation)

        # --- CTCF relative BED + overlay figure ---
        positions: list[int] = []
        orientations: list[int] = []
        if do_ctcf_overlay and ctcf_bed_path and os.path.exists(ctcf_bed_path):
            try:
                _, positions, orientations = ctcf_mod.emit_ctcf_bed_and_figure(
                    sim_map,
                    bed_path=ctcf_bed_path,
                    out_bed=os.path.join(common_dir, f"{file_prefix}_ctcf_sites_relative.bed"),
                    out_fig=os.path.join(common_dir, f"{file_prefix}_contact_map_with_ctcf.png"),
                    chrom=CHROM,
                    region_start_bp=REGION_START,
                    region_end_bp=REGION_END,
                    resolution_bp=RESOLUTION,
                    title=f"{base_name}: {cell_type or ''} contact map + CTCF",
                    elements_bed=elements_bed_path,
                    elements_label=elements_label,
                )
            except Exception as e:
                logger.warning(f"    CTCF overlay failed: {e}")

        # --- Absolute loop quantification (APA + per-loop strength) ---
        if do_apa and positions:
            try:
                loops = abs_mod.enumerate_convergent_ctcf_loops(
                    positions, orientations,
                    min_size_bins=20, max_size_bins=min(sim_map.shape[0], 1000),
                )
                if loops:
                    pileup, n_used = abs_mod.apa_pileup(
                        sim_map, loops, pad=10, normalise="obs_over_exp_diag",
                    )
                    abs_mod.plot_apa(
                        pileup,
                        os.path.join(common_dir, f"{file_prefix}_apa_convergent.png"),
                        title=f"{base_name}: APA (n={n_used} conv. CTCF pairs)",
                        n_loops=n_used,
                    )
                    quants = abs_mod.batch_absolute_quant(sim_map, loops)
                    with open(os.path.join(common_dir, f"{file_prefix}_apa_loops_quant.json"), "w") as f:
                        json.dump({"loops": quants, "n_pairs": n_used}, f, indent=2)
            except Exception as e:
                logger.warning(f"    APA/AbLE failed: {e}")

        # --- Experimental comparison ---
        metrics: dict = {}
        if hic_path and os.path.exists(hic_path):
            exp_map = np.load(hic_path)
            metrics = compare_contact_maps(sim_map, exp_map)
            with open(os.path.join(common_dir, f"{file_prefix}_metrics.json"), "w") as f:
                json.dump(metrics, f, indent=2)
            logger.info(f"    {base_name}: SCC={metrics['stratum_adjusted_corr']:.4f}  "
                        f"Pearson={metrics['overall_pearson']:.4f}")
        elif mcool_path:
            exp_map = exp_mod.load_region_matrix(
                mcool_path, CHROM, REGION_START, REGION_END, RESOLUTION,
            )
            if exp_map is not None:
                m = exp_mod.compute_sim_exp_metrics(sim_map, exp_map)
                exp_mod.save_metrics_json(
                    m, os.path.join(common_dir, f"{file_prefix}_metrics.json"),
                )
                exp_mod.plot_sim_vs_exp_map(
                    sim_map, exp_map,
                    os.path.join(common_dir, f"{file_prefix}_sim_vs_exp_map.png"),
                    sim_label=f"sim: {base_name}",
                    exp_label=f"exp: {os.path.basename(mcool_path)}",
                )
                logger.info(f"    {base_name}: (mcool) SCC={m['stratum_adjusted_corr']:.4f}")

        # --- Two-point MSD + polymer dynamics + calibration ---
        # In resume mode, skip the MSD/polymer-dynamics block entirely if the
        # canonical cached outputs (msd_overlay.png + rg_timecourse.npz) are
        # already present. Both are written at the very end of the block, so
        # their presence is a reliable "done" marker.
        msd_cached = (os.path.exists(os.path.join(common_dir, f"{file_prefix}_msd_overlay.png"))
                      and os.path.exists(os.path.join(common_dir, f"{file_prefix}_rg_timecourse.npz")))

        if reuse_heavy and msd_cached and (do_msd or do_polymer_dynamics):
            logger.info(f"    {base_name}: [resume] reusing cached MSD / "
                        f"polymer-dynamics outputs")
        elif (do_msd or do_polymer_dynamics):
            # If we took the cached-contact-map shortcut above, conformations
            # were never loaded. Load them now — MSD / dynamics need them.
            if conformations is None:
                try:
                    conformations = load_conformations_h5(sim_dir)
                except FileNotFoundError:
                    conformations = None

            if conformations is None:
                logger.info(f"    {base_name}: no conformations available for "
                            f"MSD / polymer-dynamics; skipping.")
            else:
                try:
                    # Try to grab experimental P(s) for the "hic" calibration anchor
                    sim_ps_curve = (distances, ps)
                    expt_ps_curve = None
                    if calibrate_with == "hic" and mcool_path and os.path.exists(mcool_path):
                        try:
                            e_d, e_p = ps_mod.compute_ps_from_cooler(
                                mcool_path, region=f"{CHROM}:{REGION_START}-{REGION_END}",
                                resolution=RESOLUTION,
                            )
                            expt_ps_curve = (e_d, e_p)
                        except Exception as e:
                            logger.info(f"    [calib hic] could not compute expt P(s): {e}")

                    _run_msd_and_dynamics(
                        conformations, cell_type=cell_type or CTCF_TYPES[0],
                        common_dir=common_dir, file_prefix=file_prefix,
                        display_name=base_name,
                        do_msd=do_msd,
                        do_polymer_dynamics=do_polymer_dynamics,
                        calibrate_with=calibrate_with,
                        expt_msd_csv=expt_msd_csv,
                        expt_ps_curve=expt_ps_curve,
                        sim_ps_curve=sim_ps_curve,
                    )
                except Exception as e:
                    logger.warning(f"    MSD / polymer-dynamics block failed: {e}")

        if cleanup_legacy_cache and os.path.isdir(legacy_dir):
            try:
                shutil.rmtree(legacy_dir)
                logger.info(f"  [legacy cache] removed {legacy_dir}")
            except Exception as e:  # noqa: BLE001
                logger.warning(f"  [legacy cache] failed to remove {legacy_dir}: {e}")

        logger.info(f"  Done: {base_name} → {common_dir}/{file_prefix}_*")
        return True

    except Exception:
        logger.error(f"  FAILED: {base_name}\n{traceback.format_exc()}")
        return False


# ---------------------------------------------------------------------------
# Pooled replicate analysis (Hansen lab style)
# ---------------------------------------------------------------------------

def pool_replicates(results_dir: str, common_dir: str, hic_dir: str | None,
                    n_jobs: int,
                    ctcf_bed_overrides: dict[str, str] | None = None,
                    mcool_overrides: dict[str, str] | None = None,
                    elements_bed_overrides: dict[str, str] | None = None,
                    elements_label: str = "enhancers/promoters",
                    do_apa: bool = True,
                    do_ctcf_overlay: bool = True,
                    do_msd: bool = True,
                    do_polymer_dynamics: bool = True,
                    calibrate_with: str = "hic",
                    expt_msd_overrides: dict[str, str] | None = None,
                    repo_root: str | None = None,
                    reuse_heavy: bool = False,
                    merged_dirs: list[str] | None = None):
    """
    Pool conformations across all replicates for each condition, then compute
    a single contact map from the combined ensemble and run the full modular
    analysis on that pooled map (P(s), insulation, CTCF overlay, APA, exp
    comparison).

    Also writes ``ps_overlay_all_conditions.png`` summarising every pooled
    condition on one log-log axis.

    ``merged_dirs`` limits pooling to those merged directories (used by
    ``--condition``); by default every merged dir in ``results_dir`` is pooled.
    """
    import numpy as np
    from contact_maps import (
        load_conformations_h5, load_lef_contact_map,
        extract_tiles_and_average,
        compute_insulation_score,
        compare_contact_maps,
    )
    import ps_curve as ps_mod
    import absolute_quant as abs_mod
    import ctcf_plotting as ctcf_mod
    import experimental_compare as exp_mod

    from configs.parameters import RESOLUTION, CHROM, REGION_START, REGION_END
    repo_root = repo_root or _REPO_ROOT

    pooled_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    pooled_maps: dict[str, np.ndarray] = {}

    # Group directories by condition (strip _rep<N>)
    all_dirs = merged_dirs if merged_dirs is not None else discover_merged_dirs(results_dir)
    condition_groups = defaultdict(list)
    for d in all_dirs:
        base = os.path.basename(d)
        condition, rep = parse_condition_rep(base)
        if rep >= 0:
            condition_groups[condition].append((rep, d))

    if not condition_groups:
        logger.info("No replicate groups found to pool.")
        return

    for condition in sorted(condition_groups):
        reps = sorted(condition_groups[condition])
        rep_dirs = [d for _, d in reps]
        rep_nums = [r for r, _ in reps]
        logger.info(f"  Pooling {condition}: {len(reps)} replicates "
                    f"(reps {rep_nums})")

        try:
            rep_streams: list = []
            all_conformations = None   # assigned from rep_streams after the load loop
            total_blocks = 0
            has_lef_only = False
            sim_map = None
            prefix = None

            # --- Resume-mode fast path ---
            # If a pooled contact map with the same {condition}_*blk_pooled stem
            # already exists in common_dir, reuse it instead of re-extracting
            # contacts from the full conformation ensemble.
            if reuse_heavy:
                cached = sorted(glob.glob(os.path.join(
                    common_dir, f"{condition}_*blk_pooled_contact_map.npy")))
                if cached:
                    cached_map_path = cached[-1]
                    sim_map = np.load(cached_map_path)
                    prefix = os.path.basename(cached_map_path).rsplit(
                        "_contact_map.npy", 1)[0]
                    # Parse total_blocks out of the cached filename so we keep
                    # consistent stems for downstream outputs.
                    m = re.search(r"_(\d+)blk_pooled$", prefix)
                    total_blocks = int(m.group(1)) if m else 0
                    logger.info(f"    {condition}: [resume] loaded cached "
                                f"pooled contact map ({prefix})")

            if sim_map is None:
                # --- Open conformations from ALL replicates as streams ---
                # ConformationStream objects hold no frame data in memory;
                # each pass re-reads from disk. This replaces the old
                # `all_conformations.extend(confs)` pattern which materialised
                # every replicate's trajectory into one list (~500 GB for 3
                # reps of 100k-frame 70k-monomer sims → OOM in Phase 2).
                for rep_num, rep_dir in reps:
                    try:
                        confs = load_conformations_h5(rep_dir)
                        if confs is not None:
                            n = len(confs) if hasattr(confs, '__len__') else 0
                            total_blocks += n
                            rep_streams.append(confs)
                            logger.info(f"    rep{rep_num}: {n} conformations (streaming)")
                        else:
                            # LEF-only: load contact map directly
                            has_lef_only = True
                            logger.info(f"    rep{rep_num}: LEF contact map")
                    except FileNotFoundError:
                        logger.warning(f"    rep{rep_num}: no data found, skipping")
                        continue

                if not rep_streams and not has_lef_only:
                    logger.warning(f"  {condition}: no conformations found, skipping pooling")
                    continue

                # Chain the per-replicate streams into one re-iterable pool.
                if rep_streams:
                    from analysis.contact_maps import ConformationStream
                    all_conformations = (rep_streams[0] if len(rep_streams) == 1
                                         else ConformationStream.chain(*rep_streams))

                # --- Compute pooled contact map ---
                if all_conformations is not None:
                    from configs.parameters import SIM_RUN
                    contact_radius = SIM_RUN["contact_radius"]
                    logger.info(f"    {condition}: {total_blocks} total conformations "
                                f"pooled (streaming from {len(rep_streams)} reps), {n_jobs} workers")
                    sim_map = extract_tiles_and_average(
                        all_conformations, contact_radius, n_jobs
                    )
                else:
                    # LEF-only: average the per-rep contact maps
                    lef_maps = []
                    for _, rep_dir in reps:
                        lef_map = load_lef_contact_map(rep_dir)
                        if lef_map is not None:
                            lef_maps.append(lef_map)
                            total_blocks += 1
                    if not lef_maps:
                        logger.warning(f"  {condition}: no LEF maps found, skipping")
                        continue
                    sim_map = np.mean(lef_maps, axis=0)

                # --- Save pooled outputs ---
                prefix = f"{condition}_{total_blocks}blk_pooled"

                np.save(os.path.join(common_dir, f"{prefix}_contact_map.npy"), sim_map)

            # At this point sim_map and prefix are guaranteed set (either from
            # the cache fast-path or the compute path above). P(s) and
            # insulation are cheap — always recompute so they stay in sync
            # with the contact map that will be used by the downstream steps.
            distances, ps = ps_mod.compute_ps_curve(sim_map, RESOLUTION)
            np.savez(os.path.join(common_dir, f"{prefix}_ps_curve.npz"),
                     distances=distances, ps=ps)
            pooled_curves[condition] = (distances, ps)
            pooled_maps[condition] = sim_map

            insulation = compute_insulation_score(sim_map)
            np.save(os.path.join(common_dir, f"{prefix}_insulation.npy"), insulation)

            # --- P(s) metrics ---
            ps_payload = ps_mod.save_ps_json(
                os.path.join(common_dir, f"{prefix}_ps_metrics.json"),
                distances, ps, RESOLUTION,
            )
            exp_alpha = ps_payload["metrics"].get("ps_exponent_2_30kb")
            if exp_alpha is not None:
                logger.info(f"    {condition} pooled: P(s) exponent α = {exp_alpha:.3f}")

            # --- CTCF relative BED + overlay on the pooled map ---
            cell_type, bed_path = resolve_ctcf_bed(
                condition, ctcf_bed_overrides or {}, repo_root,
            )
            elements_path = resolve_elements_bed(
                condition, elements_bed_overrides or {},
            )
            positions: list[int] = []
            orientations: list[int] = []
            if do_ctcf_overlay and bed_path and os.path.exists(bed_path):
                try:
                    _, positions, orientations = ctcf_mod.emit_ctcf_bed_and_figure(
                        sim_map,
                        bed_path=bed_path,
                        out_bed=os.path.join(common_dir,
                                             f"{prefix}_ctcf_sites_relative.bed"),
                        out_fig=os.path.join(common_dir,
                                             f"{prefix}_contact_map_with_ctcf.png"),
                        chrom=CHROM,
                        region_start_bp=REGION_START,
                        region_end_bp=REGION_END,
                        resolution_bp=RESOLUTION,
                        title=f"{condition} pooled: {cell_type}",
                        elements_bed=elements_path,
                        elements_label=elements_label,
                    )
                except Exception as e:
                    logger.warning(f"    pooled CTCF overlay failed: {e}")

            # --- Pooled APA + AbLE ---
            if do_apa and positions:
                try:
                    loops = abs_mod.enumerate_convergent_ctcf_loops(
                        positions, orientations,
                        min_size_bins=20,
                        max_size_bins=min(sim_map.shape[0], 1000),
                    )
                    if loops:
                        pileup, n_used = abs_mod.apa_pileup(
                            sim_map, loops, pad=10,
                            normalise="obs_over_exp_diag",
                        )
                        abs_mod.plot_apa(
                            pileup,
                            os.path.join(common_dir,
                                         f"{prefix}_apa_convergent.png"),
                            title=f"{condition} pooled: APA (n={n_used})",
                            n_loops=n_used,
                        )
                        quants = abs_mod.batch_absolute_quant(sim_map, loops)
                        with open(os.path.join(common_dir,
                                               f"{prefix}_apa_loops_quant.json"),
                                  "w") as f:
                            json.dump({"loops": quants, "n_pairs": n_used},
                                      f, indent=2)
                except Exception as e:
                    logger.warning(f"    pooled APA/AbLE failed: {e}")

            # --- Experimental comparison: .npy first, then mcool ---
            hic_path = find_hic_path(condition, hic_dir)
            mcool_path = resolve_mcool(condition, hic_dir, mcool_overrides or {})

            if hic_path and os.path.exists(hic_path):
                exp_map = np.load(hic_path)
                metrics = compare_contact_maps(sim_map, exp_map)
                with open(os.path.join(common_dir, f"{prefix}_metrics.json"), "w") as f:
                    json.dump(metrics, f, indent=2)
                logger.info(f"    {condition} pooled: SCC={metrics['stratum_adjusted_corr']:.4f}  "
                            f"Pearson={metrics['overall_pearson']:.4f}")
            elif mcool_path:
                exp_map = exp_mod.load_region_matrix(
                    mcool_path, CHROM, REGION_START, REGION_END, RESOLUTION,
                )
                if exp_map is not None:
                    m = exp_mod.compute_sim_exp_metrics(sim_map, exp_map)
                    exp_mod.save_metrics_json(
                        m, os.path.join(common_dir, f"{prefix}_metrics.json"),
                    )
                    exp_mod.plot_sim_vs_exp_map(
                        sim_map, exp_map,
                        os.path.join(common_dir, f"{prefix}_sim_vs_exp_map.png"),
                        sim_label=f"{condition} pooled",
                        exp_label=os.path.basename(mcool_path),
                    )
                    logger.info(f"    {condition} pooled: (mcool) "
                                f"SCC={m['stratum_adjusted_corr']:.4f}")
            else:
                logger.info(f"    {condition} pooled: no experimental reference")

            # --- Pooled MSD + polymer dynamics + calibration ---
            # Resume mode: skip if the canonical cached outputs are already
            # present in common_dir under the {prefix}_ stem.
            pooled_msd_cached = (
                os.path.exists(os.path.join(common_dir, f"{prefix}_msd_overlay.png"))
                and os.path.exists(os.path.join(common_dir, f"{prefix}_rg_timecourse.npz"))
            )
            if reuse_heavy and pooled_msd_cached and (do_msd or do_polymer_dynamics):
                logger.info(f"    {condition} pooled: [resume] reusing cached "
                            f"MSD / polymer-dynamics outputs")
            elif (do_msd or do_polymer_dynamics) and all_conformations:
                try:
                    expt_ps_curve = None
                    if calibrate_with == "hic" and mcool_path and os.path.exists(mcool_path):
                        try:
                            e_d, e_p = ps_mod.compute_ps_from_cooler(
                                mcool_path,
                                region=f"{CHROM}:{REGION_START}-{REGION_END}",
                                resolution=RESOLUTION,
                            )
                            expt_ps_curve = (e_d, e_p)
                        except Exception as e:
                            logger.info(f"    [calib hic] expt P(s) unavailable: {e}")

                    expt_msd_csv = (expt_msd_overrides or {}).get(cell_type)

                    _run_msd_and_dynamics(
                        all_conformations, cell_type=cell_type,
                        common_dir=common_dir, file_prefix=prefix,
                        display_name=f"{condition} pooled",
                        do_msd=do_msd,
                        do_polymer_dynamics=do_polymer_dynamics,
                        calibrate_with=calibrate_with,
                        expt_msd_csv=expt_msd_csv,
                        expt_ps_curve=expt_ps_curve,
                        sim_ps_curve=(distances, ps),
                    )
                except Exception as e:
                    logger.warning(f"    pooled MSD / dynamics failed: {e}")

            logger.info(f"  Done pooling: {condition} → {prefix}_*")

        except Exception:
            logger.error(f"  FAILED pooling: {condition}\n{traceback.format_exc()}")

    # --- Overlay plot of all pooled P(s) curves ---
    if pooled_curves:
        overlay_path = os.path.join(common_dir, "ps_overlay_all_conditions.png")
        try:
            ps_mod.plot_ps_overlay(
                pooled_curves, overlay_path,
                title="P(s): all pooled conditions",
            )
        except Exception as e:
            logger.warning(f"  P(s) overlay failed: {e}")

        # --- Comparative first-derivative analysis of pooled P(s) ---
        # d log P / d log s is the local contact-scaling exponent and is
        # often a cleaner between-condition discriminator than the full
        # P(s) curve. Experimental P(s) curves are also pulled from the
        # provided mcool files (once per cell type) so the derivative
        # plot and table directly compare sim vs. in-vivo. See README
        # section 6.
        try:
            deriv_png = os.path.join(
                common_dir, "ps_derivative_all_conditions.png")
            deriv_csv = os.path.join(
                common_dir, "ps_derivative_table.csv")
            deriv_json = os.path.join(
                common_dir, "ps_derivative_table.json")

            # Pull experimental P(s) from each Hi-C map once.
            exp_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
            for ct_key in MCOOL_FILES:
                mc_path = find_mcool(ct_key, hic_dir, mcool_overrides or {})
                if mc_path:
                    try:
                        e_d, e_p = ps_mod.compute_ps_from_cooler(
                            mc_path,
                            region=f"{CHROM}:{REGION_START}-{REGION_END}",
                            resolution=RESOLUTION,
                        )
                        exp_label = f"{ct_key}_experimental (exp)"
                        exp_curves[exp_label] = (e_d, e_p)
                        logger.info(
                            f"    derivative: loaded experimental P(s) "
                            f"for {ct_key} from {os.path.basename(mc_path)}"
                        )
                    except Exception as e:
                        logger.info(
                            f"    derivative: could not load experimental "
                            f"P(s) for {ct_key} ({e}); skipping exp row."
                        )

            combined_curves = {**pooled_curves, **exp_curves}

            ps_mod.plot_ps_derivative_overlay(
                combined_curves, deriv_png,
                title="Local slope of P(s): sim + experimental",
            )
            table = ps_mod.summarize_ps_derivatives(combined_curves)
            ps_mod.save_ps_derivative_table(
                table, csv_path=deriv_csv, json_path=deriv_json,
            )
        except Exception as e:
            logger.warning(f"  P(s) derivative analysis failed: {e}")

    # --- AbLE on the conserved MSD CTCF pairs, across all pooled conditions ---
    # See README section 7. This uses the same pair list picked by
    # configs/parameters.get_msd_pairs() so MSD changes and loop-strength
    # changes are always reported on identical anchors.
    if pooled_maps:
        try:
            from configs.parameters import get_msd_pairs
            # With require_conserved=True the pairs are identical across
            # ctcf_types, so using the first type is enough.
            msd_pairs, _ = get_msd_pairs(CTCF_TYPES[0])
        except Exception as e:
            logger.warning(f"  AbLE: could not resolve MSD pairs ({e}); "
                           f"skipping AbLE-on-conserved-pairs table.")
            msd_pairs = []

        if msd_pairs:
            try:
                able_res = abs_mod.able_pairs_across_conditions(
                    pooled_maps, msd_pairs,
                )
                abs_mod.save_able_table(
                    able_res,
                    csv_path=os.path.join(
                        common_dir, "able_conserved_pairs_table.csv"),
                    json_path=os.path.join(
                        common_dir, "able_conserved_pairs_table.json"),
                )
                abs_mod.plot_able_heatmap(
                    able_res,
                    out_path=os.path.join(
                        common_dir, "able_conserved_pairs_heatmap.png"),
                    title=("AbLE strength on conserved MSD pairs "
                           "(observed minus local P(s))"),
                )
            except Exception as e:
                logger.warning(
                    f"  AbLE-on-conserved-pairs analysis failed: {e}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run contact map analysis on all merged simulation directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--results-dir", type=str, required=True,
                        help="Root directory containing merged simulation dirs")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Common output folder for all analysis results. "
                             "Default: results/analysis/ (sibling of --results-dir)")
    parser.add_argument("--condition", type=str, default=None,
                        help="Analyse only this specific merged dir name "
                             "(a condition name from the config). Default: all.")
    parser.add_argument("--hic-dir", type=str, default=None,
                        help="Directory with experimental Hi-C named after the "
                             "config's compare_hic types (<type>.cool or <type>.npy). "
                             "Falls back to data/mcool/<type>.cool.")
    parser.add_argument("--n-jobs", type=int, default=4,
                        help="Parallel workers for contact detection within each dir "
                             "(default 4). Each directory analysis uses this many cores.")
    parser.add_argument("--no-parallel", action="store_true",
                        help="Analyse directories sequentially instead of in parallel. "
                             "Useful for debugging or when RAM is tight.")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip directories that already have "
                             "analysis/sim_contact_map.npy (COARSE: skips the "
                             "entire per-replicate analysis, including APA / "
                             "CTCF overlay / MSD). For a resume pass that keeps "
                             "upstream cached but redoes downstream, use "
                             "--reuse-heavy instead.")
    parser.add_argument("--reuse-heavy", action="store_true",
                        help="Resume mode: when sim_contact_map.npy (per "
                             "replicate) or {condition}_*_pooled_contact_map.npy "
                             "(pooled) already exists, load it from disk "
                             "instead of re-extracting contacts from "
                             "conformations. Likewise reuse MSD / polymer-"
                             "dynamics outputs when their canonical .png / "
                             ".npz markers are present. P(s), insulation, APA, "
                             "CTCF overlay, experimental comparison, summary "
                             "tables and Phase 3/4 are ALWAYS re-run so new "
                             "CTCF BEDs or settings take effect. Intended for "
                             "the second run after changing only CTCF inputs "
                             "or the downstream cosmetic/quantification steps.")
    parser.add_argument("--cleanup-legacy-cache", action="store_true",
                        help="After each successful per-rep run, delete the "
                             "legacy <sim_dir>/analysis/ directory (which is "
                             "no longer used by the flat-folder layout).")
    parser.add_argument("--no-pool", action="store_true",
                        help="Skip the pooled-replicate step. Only run per-rep analysis.")
    parser.add_argument("--no-apa", action="store_true",
                        help="Skip the APA / absolute loop quantification step.")
    parser.add_argument("--no-ctcf-overlay", action="store_true",
                        help="Skip the CTCF relative BED + aligned figure step.")
    parser.add_argument("--mcool", action="append", metavar="HIC_TYPE=PATH",
                        help="Override the experimental mcool for a compare_hic "
                             "type from the config. Repeatable.")
    parser.add_argument("--ctcf-bed", action="append", metavar="CTCF_TYPE=PATH",
                        help="Override the oriented CTCF BED for a ctcf_type. Repeatable.")
    parser.add_argument("--elements-bed", action="append", metavar="CTCF_TYPE=PATH",
                        help="Optional BED of non-CTCF sticky elements "
                             "(enhancers/promoters/boundaries) for conditions "
                             "with that ctcf_type. Drawn as '■' markers on the "
                             "1D tracks above and left of the contact map. Repeatable.")
    parser.add_argument("--elements-label", type=str, default="enhancers/promoters",
                        help="Legend label for the sticky-element overlay "
                             "(default: 'enhancers/promoters').")
    # --- two-point MSD + polymer dynamics + calibration ---
    parser.add_argument("--no-msd", action="store_true",
                        help="Skip the two-point MSD analysis (fastest exit "
                             "from the new block, keeps polymer-dynamics on "
                             "unless you also pass --no-polymer-dynamics).")
    parser.add_argument("--no-polymer-dynamics", action="store_true",
                        help="Skip R_g / dwell-time / looped-fraction / "
                             "pair-distance analyses.")
    parser.add_argument("--calibrate-with", choices=("hic", "msd", "none"),
                        default="hic",
                        help="Physical-units anchor strategy. "
                             "'hic': match sim P(s) to experimental Hi-C P(s) "
                             "at a reference separation (default; uses the "
                             "mcool files already provided above). "
                             "'msd': Fbn2-style MSD match, needs --expt-msd. "
                             "'none': skip calibration.")
    parser.add_argument("--expt-msd", action="append", metavar="CTCF_TYPE=PATH",
                        help="CSV with columns (dt_s, msd_um2) giving the "
                             "experimental two-point MSD at the simulated locus "
                             "for that ctcf_type; used by --calibrate-with msd. "
                             "Repeatable.")
    # --- between-condition MSD statistics (Phase 3) ---
    parser.add_argument("--no-msd-stats", action="store_true",
                        help="Skip the between-condition MSD statistics step "
                             "(Phase 3). Only relevant when multiple "
                             "conditions are analysed together.")
    parser.add_argument("--n-boot", type=int, default=10_000,
                        help="Number of bootstrap resamples AND permutation "
                             "resamples for MSD statistics (default: 10000).")
    parser.add_argument("--stats-seed", type=int, default=0,
                        help="Random seed for bootstrap and permutation tests "
                             "(default: 0, reproducible).")
    parser.add_argument("--min-reps-for-stats", type=int, default=2,
                        help="Minimum number of successful per-replicate fits "
                             "required to include a condition in the "
                             "between-condition comparisons (default: 2; "
                             "raise to 3 for production).")
    # --- Phase 4: single pivoted summary table across all conditions ---
    parser.add_argument("--no-summary-table", action="store_true",
                        help="Skip the final metric-by-condition summary "
                             "table (Phase 4). Normally this step collects "
                             "every scalar produced by the previous phases "
                             "into a single CSV/JSON/XLSX for at-a-glance "
                             "comparison.")
    parser.add_argument("--no-stats-in-summary", action="store_true",
                        help="In Phase 4, skip the significance-testing "
                             "companion CSVs (summary_stats_all_pairs and "
                             "summary_stats_vs_reference). The pivoted "
                             "summary itself is still written.")
    parser.add_argument("--summary-reference-condition", type=str, default=None,
                        help="Condition to use as the reference column in "
                             "summary_stats_vs_reference.csv. Defaults to "
                             "the first condition discovered (alphabetical "
                             "within the pooled set).")
    args = parser.parse_args()
    mcool_overrides    = parse_type_overrides(args.mcool, "--mcool")
    ctcf_bed_overrides = parse_type_overrides(args.ctcf_bed, "--ctcf-bed")
    elements_overrides = parse_type_overrides(args.elements_bed, "--elements-bed")
    expt_msd_overrides = parse_type_overrides(args.expt_msd, "--expt-msd")

    if not os.path.isdir(args.results_dir):
        logger.error(f"Results directory not found: {args.results_dir}")
        sys.exit(1)

    # Common output folder
    if args.output_dir:
        common_dir = args.output_dir
    else:
        common_dir = os.path.join(os.path.dirname(args.results_dir.rstrip("/")),
                                  "analysis")
    os.makedirs(common_dir, exist_ok=True)
    logger.info(f"Common output folder: {common_dir}")

    # Discover directories
    if args.condition:
        from configs.parameters import get_condition
        try:
            cond = get_condition(args.condition)
        except ValueError:
            logger.error(f"Unknown condition '{args.condition}'. "
                         f"Available: {[c['name'] for c in SIMULATION_CONDITIONS]}")
            sys.exit(1)
        params_name = cond["params"]["name"]
        ctcf_type   = cond["ctcf_type"]
        prefix = f"merged_{params_name}_ctcf-{ctcf_type}_rep"
        all_merged = discover_merged_dirs(args.results_dir)
        dirs_to_analyse = [d for d in all_merged
                           if os.path.basename(d).startswith(prefix)]
        if not dirs_to_analyse:
            logger.error(
                f"No merged dirs found for condition '{args.condition}' "
                f"(expected prefix '{prefix}' in {args.results_dir})"
            )
            sys.exit(1)
    else:
        dirs_to_analyse = discover_merged_dirs(args.results_dir)
        if not dirs_to_analyse:
            logger.info("No merged directories found to analyse.")
            sys.exit(0)

    logger.info(f"Found {len(dirs_to_analyse)} merged dir(s) to analyse:")
    for d in dirs_to_analyse:
        logger.info(f"  {os.path.basename(d)}")

    # Build work items (now includes common_dir)
    work_items = []
    for sim_dir in dirs_to_analyse:
        base_name = os.path.basename(sim_dir)
        hic_path  = find_hic_path(base_name, args.hic_dir)
        cell_type, ctcf_bed = resolve_ctcf_bed(
            base_name, ctcf_bed_overrides, _REPO_ROOT,
        )
        mcool_path = resolve_mcool(base_name, args.hic_dir, mcool_overrides)
        elements_bed = resolve_elements_bed(base_name, elements_overrides)
        expt_msd_csv = expt_msd_overrides.get(cell_type)
        work_items.append(dict(
            sim_dir=sim_dir,
            hic_path=hic_path,
            n_jobs=args.n_jobs,
            skip_existing=args.skip_existing,
            common_dir=common_dir,
            ctcf_bed_path=ctcf_bed,
            cell_type=cell_type,
            mcool_path=mcool_path,
            do_apa=not args.no_apa,
            do_ctcf_overlay=not args.no_ctcf_overlay,
            elements_bed_path=elements_bed,
            elements_label=args.elements_label,
            do_msd=not args.no_msd,
            do_polymer_dynamics=not args.no_polymer_dynamics,
            calibrate_with=args.calibrate_with,
            expt_msd_csv=expt_msd_csv,
            reuse_heavy=args.reuse_heavy,
            cleanup_legacy_cache=args.cleanup_legacy_cache,
        ))

    # ── Phase 1: per-replicate analysis ──────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("PHASE 1: Per-replicate analysis")
    logger.info("=" * 60)

    n_dir_workers = max(1, _allocated_cpu_count() // args.n_jobs)
    n_dir_workers = min(n_dir_workers, len(work_items))

    if not args.no_parallel and n_dir_workers > 1:
        logger.info(f"Running {len(work_items)} analyses in parallel "
                    f"({n_dir_workers} dirs × {args.n_jobs} workers each "
                    f"= {n_dir_workers * args.n_jobs} cores total)")
        results_map = {}
        with ThreadPoolExecutor(max_workers=n_dir_workers) as executor:
            futures = {
                executor.submit(run_analysis_one, **item):
                    os.path.basename(item["sim_dir"])
                for item in work_items
            }
            for future in as_completed(futures):
                name = futures[future]
                try:
                    results_map[name] = future.result()
                except Exception:
                    logger.error(f"  FAILED (uncaught): {name}\n{traceback.format_exc()}")
                    results_map[name] = False
        results = [results_map[os.path.basename(item["sim_dir"])]
                   for item in work_items]
    else:
        if args.no_parallel:
            logger.info(f"Running {len(work_items)} analyses sequentially (--no-parallel)")
        else:
            logger.info(f"Running {len(work_items)} analyses sequentially "
                        f"(single CPU available per dir)")
        results = [run_analysis_one(**item) for item in work_items]

    n_ok   = sum(results)
    n_fail = len(results) - n_ok
    logger.info(f"Phase 1 complete: {n_ok} succeeded, {n_fail} failed.")

    # ── Phase 2: pool replicates per condition ───────────────────────────
    if not args.no_pool:
        logger.info("")
        logger.info("=" * 60)
        logger.info("PHASE 2: Pooling replicates per condition (Hansen lab style)")
        logger.info("=" * 60)
        pool_replicates(
            args.results_dir, common_dir, args.hic_dir, args.n_jobs,
            ctcf_bed_overrides=ctcf_bed_overrides,
            mcool_overrides=mcool_overrides,
            elements_bed_overrides=elements_overrides,
            elements_label=args.elements_label,
            do_apa=not args.no_apa,
            do_ctcf_overlay=not args.no_ctcf_overlay,
            repo_root=_REPO_ROOT,
            do_msd=not args.no_msd,
            do_polymer_dynamics=not args.no_polymer_dynamics,
            calibrate_with=args.calibrate_with,
            expt_msd_overrides=expt_msd_overrides,
            reuse_heavy=args.reuse_heavy,
            merged_dirs=dirs_to_analyse if args.condition else None,
        )
        logger.info("Phase 2 complete.")
    else:
        logger.info("Skipping Phase 2 (--no-pool)")

    # ── Phase 3: between-condition MSD statistics ───────────────────────
    # Aggregates per-replicate alpha and K_alpha values across conditions
    # and runs pairwise comparisons (Welch's t, Mann-Whitney, permutation,
    # bootstrap CI, Cohen's d). Honest caveats on n_rep = 3 power are
    # embedded in the output JSON.
    if (not args.no_msd) and (not args.no_msd_stats) and (not args.condition):
        logger.info("")
        logger.info("=" * 60)
        logger.info("PHASE 3: Between-condition MSD statistics")
        logger.info("=" * 60)
        try:
            import msd_statistics as stats_mod
            stats_mod.run_msd_statistics(
                common_dir,
                n_boot=args.n_boot,
                seed=args.stats_seed,
                min_n_per_condition=args.min_reps_for_stats,
            )
            logger.info("Phase 3 complete.")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Phase 3 (MSD stats) failed: {e}")
    elif args.no_msd_stats:
        logger.info("Skipping Phase 3 (--no-msd-stats)")

    # ── Phase 4: metric-by-condition summary table ───────────────────────
    # Reads every scalar already written to ``common_dir`` and pivots it
    # into a single matrix (rows = metrics grouped by family, columns =
    # conditions). This is the one place to look if you want to compare
    # P(s) exponent, AbLE strength, alpha, K_alpha, SCC-vs-experiment,
    # and looped fraction side-by-side across all simulated + any
    # experimental Hi-C references pulled in by earlier phases.
    if not args.no_summary_table:
        logger.info("")
        logger.info("=" * 60)
        logger.info("PHASE 4: Metric-by-condition summary table")
        logger.info("=" * 60)
        try:
            from summarize_analysis import build_summary_table
            build_summary_table(
                common_dir,
                do_stats=not args.no_stats_in_summary,
                reference_condition=args.summary_reference_condition,
            )
            logger.info("Phase 4 complete.")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"Phase 4 (summary table) failed: {e}")
    else:
        logger.info("Skipping Phase 4 (--no-summary-table)")

    # ── Summary ─────────────────────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("ALL DONE")
    logger.info("=" * 60)
    logger.info(f"Per-replicate outputs: {n_ok} succeeded, {n_fail} failed")
    _warn_legacy_outputs(common_dir)
    logger.info(f"All outputs collected in: {common_dir}/")
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
