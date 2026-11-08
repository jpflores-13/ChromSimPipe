#!/usr/bin/env python
"""
One simulation shard: a polysim run on one GPU.

Each replicate of a condition is split into `n_shards` independent runs (each
with its own seed); merge_shards.py concatenates their saved conformations.
The model is the Hansen lab's polysim package (polymer + loop extrusion on
polychrom / OpenMM): https://github.com/ahansenlab/polysim_bootcamp

The active locus comes from $CHROMSIM_LOCUS (set by the workflow for every
job), the settings from the config named by $CHROMSIM_CONFIG.

Usage:
    CHROMSIM_LOCUS=chr7_lost python scripts/run_simulation_shard.py \
        --condition cond1_baseline --replicate 0 \
        --shard-index 0 --n-shards 2 --gpu 0 --output <locus dir>/sims

This is called by the simulateShard rule in workflows/ChromSimPipe.snakefile.
"""

import os
import sys
import ctypes
import ctypes.util
import argparse
import logging
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import (
    SIM_RUN, POLYSIM, TILING, ACTIVE_LOCUS, get_condition, get_tiled_ctcf_arrays,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def seed_all(seed):
    """Seed numpy and libc's drand48.

    polysim's LEF translocator (LEF_Dynamics.pyx) draws its stall / unbinding /
    stepping decisions from C drand48(), which it never seeds, so without this
    every shard and replicate would replay the same random sequence.
    """
    np.random.seed(seed)
    libc = ctypes.CDLL(ctypes.util.find_library("c"))
    libc.srand48(ctypes.c_long(seed))


def polysim_kwargs(params, ctcf_pos, ctcf_ori, n_frames, gpu, shard_dir):
    """polysim SimParams keyword arguments for one shard.

    CTCF orientation follows ChromSimPipe's convention (+1 stops leftward-moving
    LEF legs, -1 stops rightward-moving ones), which is polysim's ctcf_left /
    ctcf_right. Each site stalls with the condition's ctcf_capture per encounter.
    `n_frames` is the number of conformations this shard saves.
    """
    ps = POLYSIM
    blocks = n_frames * ps["saveevery"]
    if blocks % ps["blocks_per_updater"]:
        step = ps["blocks_per_updater"] // np.gcd(ps["saveevery"], ps["blocks_per_updater"])
        raise ValueError(
            f"polysim needs frames per shard (frames_per_replicate / n_shards = {n_frames}) "
            f"times saveevery ({ps['saveevery']}) to be a multiple of blocks_per_updater "
            f"({ps['blocks_per_updater']}); use a multiple of {step} frames per shard")
    cap = float(params["ctcf_capture"])
    left = {int(p): cap for p, o in zip(ctcf_pos, ctcf_ori) if o == +1}
    right = {int(p): cap for p, o in zip(ctcf_pos, ctcf_ori) if o == -1}
    return dict(
        npoly=int(TILING["chrom_size"]), nchr=1,
        density=ps["density"], confinement_k=ps["confinement_k"], repel=ps["repel"],
        bond_wiggle=ps["bond_wiggle"], angle_k=ps["angle_k"],
        life=float(params["lifetime"]), sep=float(params["separation"]),
        vlef=float(params.get("vlef", ps["vlef"])),
        lifebooststalled=float(params.get("life_boost_stalled", ps["life_boost_stalled"])),
        stall=None, ctcf_left=left, ctcf_right=right,
        smc_bond_wiggle=ps["smc_bond_wiggle"], smc_bond_dist=ps["smc_bond_dist"],
        platform=ps["platform"], gpu=str(gpu), integrator=ps["integrator"],
        dt=ps["dt"], colrate=ps["colrate"], poly_steps_per_block=ps["poly_steps_per_block"],
        numsave=int(n_frames), saveevery=ps["saveevery"], initskip=ps["initskip"],
        initsteps=ps["initsteps"], blocks_per_updater=ps["blocks_per_updater"],
        max_data_length=SIM_RUN["hdf5_frames_per_file"], outpath=shard_dir,
        save_smc_bonds=True,
    )


def run_shard(condition, replicate, shard_index, n_shards, gpu, output_dir):
    """Run one shard and write blocks_*.h5 into <output_dir>/<name>_ctcf-<type>_rep<r>_shard<i>/."""
    cond = get_condition(condition)
    params, ctcf_type = cond["params"], cond["ctcf_type"]
    n_frames = SIM_RUN["frames_per_replicate"] // n_shards
    seed = SIM_RUN["seed_base"] + replicate * 10000 + shard_index * 1000

    shard_dir = os.path.join(output_dir,
                             f"{params['name']}_ctcf-{ctcf_type}_rep{replicate}_shard{shard_index}")
    os.makedirs(shard_dir, exist_ok=True)
    with open(os.path.join(shard_dir, "params.json"), "w") as f:
        json.dump({**params, "condition": condition, "ctcf_type": ctcf_type,
                   "locus": ACTIVE_LOCUS, "replicate": replicate,
                   "shard_index": shard_index, "n_shards": n_shards,
                   "frames": n_frames, "seed": seed, "tiling": TILING,
                   "polysim": POLYSIM}, f, indent=2, default=str)

    try:
        from polysim.sim3d import SimParams, run
    except ImportError as e:
        raise SystemExit(f"polysim cannot be imported ({e}). Rebuild the condaEnvs.sim environment.")

    seed_all(seed)
    ctcf_pos, ctcf_ori = get_tiled_ctcf_arrays(ctcf_type)
    sp = SimParams(**polysim_kwargs(params, ctcf_pos, ctcf_ori, n_frames, gpu, shard_dir))
    logger.info(f"{ACTIVE_LOCUS} / {condition} rep{replicate} shard{shard_index}: "
                f"{n_frames} frames, seed {seed}, GPU {gpu}")
    run(sp, folder=shard_dir, verbose=True)
    logger.info("Shard complete")
    return shard_dir


def main():
    parser = argparse.ArgumentParser(description="Run one polysim simulation shard")
    parser.add_argument("--condition", required=True, help="Condition name from the config")
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--output", required=True, help="The locus' sims folder")
    args = parser.parse_args()
    run_shard(args.condition, args.replicate, args.shard_index, args.n_shards,
              args.gpu, args.output)


if __name__ == "__main__":
    main()
