#!/usr/bin/env python
"""
Fast-polymer-dynamics variant of ``run_analysis_all.py``.

Identical CLI, identical on-disk outputs, identical Phase 1 / Phase 2 /
Phase 3 / Phase 4 orchestration — the **only** difference is that the
polymer-dynamics block inside ``_run_msd_and_dynamics`` is replaced with
a single file-parallel pass provided by
``analysis.polymer_dynamics_parallel.compute_dynamics_batch``.

Why a new script
----------------
The current ``_run_msd_and_dynamics`` runs polymer dynamics as seven
sequential single-threaded Python loops over the pooled trajectory
(R_g + looped_fraction × n_pairs + dwell_times × n_pairs +
pair_distance_dist × n_pairs). On pooled 3-replicate conditions this
costs ~2 h per condition while 15 of 16 allocated cores idle.

This script routes that same work through the Hansen-lab file-parallel
worker pattern already used for contact-map extraction, so the
polymer-dynamics block drops from ~2 h to ~10 min per condition.

It does NOT touch ``scripts/run_analysis_all.py`` or
``analysis/polymer_dynamics.py`` — both remain available for regression
comparison.

Usage
-----
Drop-in replacement for ``run_analysis_all.py``::

    python scripts/run_analysis_all_fastdyn.py \
        --results-dir results/polychrom_3d \
        --n-jobs 16 \
        --reuse-heavy --skip-existing \
        ... (every flag accepted by run_analysis_all.py is accepted here)

From a SLURM job, point the existing ``submit_analysis_all_16cpu.sh``
script at this file by setting ``ANALYSIS_ENTRY=run_analysis_all_fastdyn``
(see ``cluster/submit_analysis_all_16cpu_fastdyn.sh`` for the matching
sbatch wrapper).

Parallelism
-----------
The batched polymer-dynamics call reads the number of workers from, in
priority order:

    DYN_N_JOBS              -> explicit override (rarely needed)
    SLURM_CPUS_PER_TASK     -> normal cluster path (matches --n-jobs)
    os.cpu_count() or 16    -> final fallback

Set ``DYN_N_JOBS`` if you want the polymer-dynamics fan-out smaller
than the MSD one (for RAM-tight nodes); leaving it unset mirrors the
cluster's allocation.

Failure mode
------------
If a trajectory has no ``block_files`` populated (legacy
``conformations.h5``, non-polychrom source), the batched function
automatically falls back to the sequential implementations in
``analysis/polymer_dynamics.py`` — so no data shape is regressed, and
the output artefacts remain the same.
"""

from __future__ import annotations

import json
import logging
import os
import sys
from typing import Optional

# --- sys.path bootstrap (mirror run_analysis_all.py) --------------------------
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
if os.path.join(_REPO_ROOT, "analysis") not in sys.path:
    sys.path.insert(0, os.path.join(_REPO_ROOT, "analysis"))

import numpy as np

# Importing run_analysis_all gives us access to its helpers, its CLI, and
# — critically — its module globals dict, which is the namespace Python
# consults when `_run_msd_and_dynamics` is referenced bare inside
# `analyze_one_replicate` / `pool_replicates`. Monkey-patching an entry
# in that dict is enough to redirect every call.
import run_analysis_all as ram  # noqa: E402

logger = logging.getLogger(__name__)


# =============================================================================
# HELPERS
# =============================================================================

def _get_dyn_njobs() -> int:
    """Resolve the worker count used for the parallel polymer-dynamics pass.

    Priority: DYN_N_JOBS > SLURM_CPUS_PER_TASK > os.cpu_count() > 16.
    """
    for var in ("DYN_N_JOBS", "SLURM_CPUS_PER_TASK"):
        v = os.environ.get(var, "").strip()
        if v:
            try:
                n = int(v)
                if n > 0:
                    return n
            except ValueError:
                pass
    n = os.cpu_count() or 16
    return max(1, n)


# =============================================================================
# REPLACEMENT _run_msd_and_dynamics
# =============================================================================
#
# Byte-for-byte faithful to ``run_analysis_all._run_msd_and_dynamics`` with
# the exception of the polymer-dynamics block, which is collapsed into one
# call to ``compute_dynamics_batch``. The MSD loop, MSD overlay plot, and
# physical-units calibration are preserved verbatim.
#
# The signature MUST match the original so that the call sites inside
# ``analyze_one_replicate`` (line ~861) and ``pool_replicates`` (line ~1233)
# keep working after the monkey-patch.
# =============================================================================

def _run_msd_and_dynamics_fastdyn(
    conformations,
    cell_type: str,
    out_dir: str,
    display_name: str = "",
    *,
    do_msd: bool = True,
    do_polymer_dynamics: bool = True,
    calibrate_with: str = "hic",
    expt_msd_csv: Optional[str] = None,
    expt_ps_curve: Optional[tuple] = None,
    sim_ps_curve: Optional[tuple] = None,
) -> dict:
    """Same contract as ``run_analysis_all._run_msd_and_dynamics`` but the
    polymer-dynamics block runs file-parallel via
    ``polymer_dynamics_parallel.compute_dynamics_batch``.
    """
    # Short-circuit on empty / LEF-only inputs — matches the original.
    if conformations is None:
        logger.info(f"    [MSD skipped: no 3D conformations for {display_name}]")
        return {"msd_curves": {}, "calibration": None,
                "pairs": [], "cell_type": cell_type}

    import msd_two_point as msd_mod
    import polymer_dynamics as pd_mod          # only for the plot helpers
    import polymer_dynamics_parallel as pdp
    import calibration as cal_mod
    from configs.parameters import (
        RESOLUTION, TILING, N_MONOMERS, SIM_RUN,
        get_msd_pairs, MSD_FIT, CALIBRATION,
    )

    tile_size = TILING["tile_size"]
    pad = TILING["padding"]
    n_tiles = TILING["n_tiles"]
    contact_radius = SIM_RUN["contact_radius"]

    os.makedirs(out_dir, exist_ok=True)
    out: dict = {
        "msd_curves": {},
        "pairs": [],
        "cell_type": cell_type,
        "calibration": None,
    }

    pairs, labels = get_msd_pairs(cell_type)
    logger.info(f"    [MSD] {cell_type}: {len(pairs)} pair(s) "
                f"{list(zip(labels, pairs))}")

    # ── Two-point MSD per pair (unchanged) ──────────────────────────────
    if do_msd:
        for (idx_a, idx_b), label in zip(pairs, labels):
            try:
                res = msd_mod.run_msd_for_pair(
                    conformations, pair=(idx_a, idx_b), label=label,
                    out_dir=out_dir,
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
            except Exception as e:                                    # noqa: BLE001
                logger.warning(f"    [MSD] pair {label} failed: {e}")

        if out["msd_curves"]:
            try:
                msd_mod.plot_msd_curves(
                    out["msd_curves"],
                    os.path.join(out_dir, "msd_overlay.png"),
                    title=f"Two-point MSD: {display_name or cell_type}",
                )
            except Exception as e:                                    # noqa: BLE001
                logger.warning(f"    [MSD] overlay failed: {e}")

    # ── Polymer-dynamics diagnostics (single file-parallel pass) ────────
    if do_polymer_dynamics:
        n_jobs = _get_dyn_njobs()
        logger.info(
            f"    [dyn] batched file-parallel pass on {display_name or cell_type} "
            f"with {n_jobs} workers (R_g + looped_fraction + dwell_times + "
            f"pair_distance_dist, {len(pairs)} pair(s))"
        )
        try:
            batch = pdp.compute_dynamics_batch(
                conformations,
                pairs=pairs,
                labels=labels,
                contact_radius=contact_radius,
                tile_size=tile_size,
                pad=pad,
                n_tiles=n_tiles,
                n_jobs=n_jobs,
                bead_range=(0, N_MONOMERS),
            )
            logger.info(
                f"    [dyn] batch source={batch['source']} "
                f"T={batch['T']} n_tiles={batch['n_tiles_used']}"
            )

            rg = batch["rg_timecourse"]

            # --- R_g artefacts (parity with the original outputs) -----
            try:
                pd_mod.plot_rg_timecourse(
                    {display_name or cell_type: rg},
                    os.path.join(out_dir, "rg_timecourse.png"),
                    title=f"Radius of gyration: {display_name or cell_type}",
                )
            except Exception as e:                                    # noqa: BLE001
                logger.warning(f"    [R_g] plot failed: {e}")

            np.savez_compressed(
                os.path.join(out_dir, "rg_timecourse.npz"),
                t=rg["t"], rg=rg["rg"],
            )
            out["rg_mean"] = float(rg["rg"].mean())
            out["rg_std"] = float(rg["rg"].std())

            # --- Per-pair artefacts ----------------------------------
            loop_fracs: list[dict] = []
            dwell_by_label: dict = {}
            pair_dist_by_label: dict = {}
            for label in labels:
                pp = batch["per_pair"][label]
                loop_fracs.append({"label": label, **pp["looped_fraction"]})
                dwell_by_label[label] = pp["dwell_times"]
                pair_dist_by_label[label] = pp["pair_distance_dist"]

            out["loop_fracs"] = loop_fracs
            if loop_fracs:
                with open(os.path.join(out_dir, "loop_fractions.json"), "w") as f:
                    json.dump(loop_fracs, f, indent=2)

            if dwell_by_label:
                try:
                    pd_mod.plot_dwell_time_hist(
                        dwell_by_label,
                        os.path.join(out_dir, "dwell_times.png"),
                        title=f"Contact dwell times: {display_name or cell_type}",
                    )
                except Exception as e:                                # noqa: BLE001
                    logger.warning(f"    [dwell] plot failed: {e}")

            if pair_dist_by_label:
                try:
                    pd_mod.plot_pair_distance_hist(
                        pair_dist_by_label,
                        os.path.join(out_dir, "pair_distance.png"),
                        title=f"Pair separation density: {display_name or cell_type}",
                    )
                except Exception as e:                                # noqa: BLE001
                    logger.warning(f"    [pair-dist] plot failed: {e}")

        except Exception as e:                                        # noqa: BLE001
            # If the batched path blows up unexpectedly, fall back to the
            # original sequential implementation so the pipeline still
            # produces artefacts for this condition. The resume pass will
            # then skip these outputs on its next iteration.
            logger.error(
                f"    [dyn] batched path failed ({e!r}); falling back to "
                f"sequential polymer_dynamics for this condition."
            )
            _sequential_polymer_dynamics_fallback(
                conformations, cell_type=cell_type, out_dir=out_dir,
                display_name=display_name,
                pairs=pairs, labels=labels,
                contact_radius=contact_radius,
                tile_size=tile_size, pad=pad, n_tiles=n_tiles,
                out=out,
            )

    # ── Physical-units calibration (unchanged) ──────────────────────────
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
                logger.info(f"    [calibration] no expt MSD CSV "
                            f"({expt_msd_csv!r}); skipping msd anchor")
            elif not out["msd_curves"]:
                logger.info("    [calibration] no sim MSD curve; "
                            "skipping msd anchor")
            else:
                first = next(iter(out["msd_curves"].values()))
                sim_lags, sim_msd = first["lags"], first["msd"]
                exp_t, exp_m = cal_mod.load_experimental_msd_csv(expt_msd_csv)
                calib = cal_mod.calibrate_from_msd(
                    sim_lags, sim_msd, exp_t, exp_m,
                )
        elif calibrate_with == "none":
            pass
        else:
            logger.warning(f"    [calibration] unknown anchor "
                           f"{calibrate_with!r}; skipping")
    except Exception as e:                                            # noqa: BLE001
        logger.warning(f"    [calibration] failed: {e}")
        calib = None

    if calib is not None:
        out["calibration"] = calib
        with open(os.path.join(out_dir, "calibration.json"), "w") as f:
            json.dump(calib.to_dict(), f, indent=2)
        logger.info(f"    [calibration] nm/monomer={calib.nm_per_monomer:.2f} "
                    f"sec/frame={calib.sec_per_frame:.4f} "
                    f"source={calib.source}")

    return out


def _sequential_polymer_dynamics_fallback(
    conformations, *, cell_type, out_dir, display_name,
    pairs, labels, contact_radius, tile_size, pad, n_tiles, out,
):
    """Emergency fallback that mirrors the original sequential block.

    Called only if ``compute_dynamics_batch`` raises. Produces the same
    on-disk artefacts so the resume pass's ``reuse_heavy`` markers
    (``msd_overlay.png`` + ``rg_timecourse.npz``) are satisfied.
    """
    import polymer_dynamics as pd_mod
    from configs.parameters import N_MONOMERS

    try:
        rg = pd_mod.compute_rg_timecourse(
            conformations, bead_range=(0, N_MONOMERS),
            tile_size=tile_size, pad=pad, n_tiles=n_tiles,
        )
        pd_mod.plot_rg_timecourse(
            {display_name or cell_type: rg},
            os.path.join(out_dir, "rg_timecourse.png"),
            title=f"Radius of gyration: {display_name or cell_type}",
        )
        np.savez_compressed(
            os.path.join(out_dir, "rg_timecourse.npz"),
            t=rg["t"], rg=rg["rg"],
        )
        out["rg_mean"] = float(rg["rg"].mean())
        out["rg_std"] = float(rg["rg"].std())
    except Exception as e:                                            # noqa: BLE001
        logger.warning(f"    [R_g fallback] failed: {e}")

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
        except Exception as e:                                        # noqa: BLE001
            logger.warning(f"    [dyn fallback] pair {label} failed: {e}")

    out["loop_fracs"] = loop_fracs
    if loop_fracs:
        with open(os.path.join(out_dir, "loop_fractions.json"), "w") as f:
            json.dump(loop_fracs, f, indent=2)
    if dwell_by_label:
        try:
            pd_mod.plot_dwell_time_hist(
                dwell_by_label,
                os.path.join(out_dir, "dwell_times.png"),
                title=f"Contact dwell times: {display_name or cell_type}",
            )
        except Exception as e:                                        # noqa: BLE001
            logger.warning(f"    [dwell fallback] plot failed: {e}")
    if pair_dist_by_label:
        try:
            pd_mod.plot_pair_distance_hist(
                pair_dist_by_label,
                os.path.join(out_dir, "pair_distance.png"),
                title=f"Pair separation density: {display_name or cell_type}",
            )
        except Exception as e:                                        # noqa: BLE001
            logger.warning(f"    [pair-dist fallback] plot failed: {e}")


# =============================================================================
# MONKEY-PATCH + ENTRY POINT
# =============================================================================

# Python resolves bare names (`_run_msd_and_dynamics(...)`) inside a
# function by consulting the defining module's globals at call time. By
# overwriting that attribute on the already-imported ``run_analysis_all``
# module we redirect every call site — both Phase 1's per-replicate path
# (``analyze_one_replicate``) and Phase 2's pooled path
# (``pool_replicates``) — without touching the original file.
ram._run_msd_and_dynamics = _run_msd_and_dynamics_fastdyn


if __name__ == "__main__":
    ram.main()
