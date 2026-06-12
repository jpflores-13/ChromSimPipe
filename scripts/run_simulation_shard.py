#!/usr/bin/env python
"""
Sharded simulation runner for multi-GPU parallelism.

Splits the total simulation blocks across N shards, each running on a
separate GPU. The LEF dynamics are seeded differently per shard to ensure
independent sampling, and each shard warms up independently.

After all shards finish, use merge_shards.py to combine conformations.

Usage:
    python run_simulation_shard.py \
        --params mESC --replicate 0 \
        --shard-index 0 --n-shards 4 --gpu 0

This is called by cluster/submit_polychrom_multigpu.sh.
"""

import os
import sys
import argparse
import logging
import json
import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import (
    POLYMER, SMC_BOND, SIM_RUN, N_MONOMERS, TILING,
    get_ctcf_arrays, get_tiled_ctcf_arrays,
)
from scripts.lef_dynamics import LEFSimulator
from scripts.run_simulation import get_param_set, resolve_condition, initialize_polymer

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def run_shard(params, replicate, shard_index, n_shards, gpu, output_dir,
              ctcf_type="mESC", shard_index_offset=0):
    """
    Run one shard of the simulation.

    The total_blocks are divided evenly across shards. Each shard gets
    a unique random seed so the conformations are independent samples
    from the same ensemble.
    """
    # --- Tiling: use the full tiled chromosome ---
    N = TILING["chrom_size"]                     # 70000
    n_tiles = TILING["n_tiles"]                  # 28
    name = params["name"]

    total_blocks = SIM_RUN["total_blocks"]
    warmup = SIM_RUN["warmup_blocks"]
    save_every = SIM_RUN["save_every"]
    md_steps = SIM_RUN["md_steps_per_block"]

    # Each shard runs (total_blocks / n_shards) production blocks
    # plus its own warmup
    blocks_per_shard = total_blocks // n_shards
    # global_shard_index accounts for any offset from previous runs,
    # ensuring unique directory names and independent random seeds
    global_shard_index = shard_index + shard_index_offset
    shard_seed = 42 + replicate * 10000 + global_shard_index * 1000

    shard_dir = os.path.join(output_dir,
                             f"{name}_ctcf-{ctcf_type}_rep{replicate}_shard{global_shard_index}")
    os.makedirs(shard_dir, exist_ok=True)

    with open(os.path.join(shard_dir, "params.json"), "w") as f:
        json.dump({
            **params, "ctcf_type": ctcf_type, "replicate": replicate,
            "shard_index": global_shard_index, "n_shards": n_shards,
            "blocks_per_shard": blocks_per_shard, "seed": shard_seed,
            "tiling": TILING,
        }, f, indent=2, default=str)

    # --- Initialize LEF simulator on the FULL tiled chromosome ---
    ctcf_pos, ctcf_ori = get_tiled_ctcf_arrays(cell_type=ctcf_type)
    n_lefs = max(1, N // params["separation"])

    lef_sim = LEFSimulator(
        N=N, n_lefs=n_lefs, lifetime=params["lifetime"],
        ctcf_positions=ctcf_pos, ctcf_orientations=ctcf_ori,
        ctcf_capture=params["ctcf_capture"], ctcf_release=params["ctcf_release"],
        rng_seed=shard_seed,
    )

    # Warm up LEFs
    logger.info(f"Shard {shard_index}: warming up LEFs ({warmup} steps)...")
    for _ in range(warmup):
        lef_sim.step()

    # --- Try polychrom, then OpenMM, then LEF-only ---
    engine_used = "lef_only"

    try:
        import polychrom
        from polychrom.simulation import Simulation
        from polychrom.forces import (
            harmonic_bonds, angle_force, polynomial_repulsive, spherical_confinement,
        )
        engine_used = "polychrom"
    except ImportError:
        try:
            import openmm
            engine_used = "openmm"
        except ImportError:
            engine_used = "lef_only"

    logger.info(f"Shard {shard_index}: engine={engine_used}, blocks={blocks_per_shard}, GPU={gpu}")

    if engine_used == "polychrom":
        _run_shard_polychrom(
            lef_sim, params, N, n_lefs, blocks_per_shard, warmup,
            save_every, md_steps, shard_dir, gpu, shard_seed
        )
    elif engine_used == "openmm":
        _run_shard_openmm(
            lef_sim, params, N, n_lefs, blocks_per_shard, warmup,
            save_every, md_steps, shard_dir, gpu, shard_seed
        )
    else:
        _run_shard_lef_only(
            lef_sim, N, blocks_per_shard, save_every, shard_dir
        )

    logger.info(f"Shard {shard_index}: complete")
    return shard_dir


def _run_shard_polychrom(lef_sim, params, N, n_lefs, n_blocks, warmup,
                         save_every, md_steps, shard_dir, gpu, seed):
    """Run polychrom-based shard.

    Uses the Hansen lab bond-update pattern from Yang JH, Brandão HB,
    Hansen AS (2023) Nat Commun 14:1913, "DNA double-strand break end
    synapsis by DNA loop extrusion" (10.1038/s41467-023-37583-w),
    `3D_PolymerSimulationCode(withLoopExtrusion).py` in
    github.com/ahansenlab/DNA_break_synapsis_models: pre-compute LEF
    states, find unique bonds, add them all to the harmonic_bonds force
    before context creation, then toggle active/inactive each block.
    """
    from polychrom.simulation import Simulation
    from polychrom.forces import (
        harmonic_bonds, angle_force, polynomial_repulsive, spherical_confinement,
    )
    from polychrom import forcekits
    from polychrom.hdf5_format import HDF5Reporter

    np.random.seed(seed)

    # --- Pre-compute all LEF bond states ---
    logger.info(f"  Pre-computing LEF bond states for {n_blocks} blocks...")
    all_bond_sets = []
    for _ in range(n_blocks):
        lef_sim.step()
        bonds = lef_sim.get_bonds()
        all_bond_sets.append([(int(l), int(r)) for l, r in bonds])

    unique_bonds = list(set(sum(all_bond_sets, [])))
    logger.info(f"  {len(unique_bonds)} unique bond pairs found")

    # --- Initialize polymer and simulation ---
    coords = initialize_polymer(N, POLYMER["density"])

    reporter = HDF5Reporter(folder=shard_dir, max_data_length=100,
                            overwrite=True, blocks_only=True)

    sim = Simulation(
        platform=POLYMER["platform"], GPU=str(gpu), N=N,
        error_tol=POLYMER["error_tol"], collision_rate=POLYMER["collision_rate"],
        integrator="variableLangevin",
        reporters=[reporter],
        PBCbox=False,
        precision="mixed",
    )
    sim.set_data(coords)

    # --- Add forces ---
    sim.add_force(spherical_confinement(
        sim, density=POLYMER["density"], k=1.0,
    ))

    sim.add_force(
        forcekits.polymer_chains(
            sim,
            chains=[(0, N, False)],  # one linear tiled chromosome
            bond_force_func=harmonic_bonds,
            bond_force_kwargs={
                "bondLength": POLYMER["bond_length"],
                "bondWiggleDistance": POLYMER["bond_wiggle_dist"],
                "override_checks": True,
            },
            angle_force_func=angle_force,
            angle_force_kwargs={
                "k": POLYMER["angle_force"],
                "override_checks": True,
            },
            nonbonded_force_func=polynomial_repulsive,
            nonbonded_force_kwargs={
                "trunc": 3.0,
            },
            except_bonds=True,
            override_checks=True,
        )
    )

    # --- Add all unique SMC bonds to the existing harmonic_bonds force ---
    kbond = sim.kbondScalingFactor / (SMC_BOND["wiggle_dist"] ** 2)
    bondDist = SMC_BOND["bond_dist"] * sim.length_scale
    activeParams = {"length": bondDist, "k": kbond}
    inactiveParams = {"length": bondDist, "k": 0}

    bondForce = sim.force_dict["harmonic_bonds"]
    first_bonds = all_bond_sets[0]
    bond_to_ind = {}
    for bond in unique_bonds:
        paramset = activeParams if (bond in first_bonds) else inactiveParams
        ind = bondForce.addBond(bond[0], bond[1], **paramset)
        bond_to_ind[bond] = ind

    # --- Create context ---
    sim.local_energy_minimization()

    # --- Main simulation loop ---
    cur_bonds = first_bonds
    for block in range(n_blocks):
        new_bonds = all_bond_sets[block]

        # Toggle changed bonds
        bonds_remove = [b for b in cur_bonds if b not in new_bonds]
        bonds_add = [b for b in new_bonds if b not in cur_bonds]
        for bond in bonds_add:
            ind = bond_to_ind[bond]
            bondForce.setBondParameters(ind, bond[0], bond[1], **activeParams)
        for bond in bonds_remove:
            ind = bond_to_ind[bond]
            bondForce.setBondParameters(ind, bond[0], bond[1], **inactiveParams)
        bondForce.updateParametersInContext(sim.context)
        cur_bonds = new_bonds

        sim.do_block(md_steps)

        if block % 500 == 0:
            logger.info(f"  block {block}/{n_blocks}")

    reporter.dump_data()
    logger.info(f"  Shard complete, {n_blocks} blocks saved")


def _run_shard_openmm(lef_sim, params, N, n_lefs, n_blocks, warmup,
                      save_every, md_steps, shard_dir, gpu, seed):
    """Run standalone OpenMM shard."""
    import openmm
    import openmm.unit as unit

    np.random.seed(seed)
    coords = initialize_polymer(N, POLYMER["density"])

    system = openmm.System()
    mass = 100.0 * unit.amu
    for _ in range(N):
        system.addParticle(mass)

    backbone = openmm.HarmonicBondForce()
    for i in range(N-1):
        backbone.addBond(i, i+1, POLYMER["bond_length"]*unit.nanometer,
                         300.0*unit.kilojoule_per_mole/unit.nanometer**2)
    system.addForce(backbone)

    angle = openmm.HarmonicAngleForce()
    for i in range(N-2):
        angle.addAngle(i, i+1, i+2, np.pi*unit.radian,
                       POLYMER["angle_force"]*unit.kilojoule_per_mole/unit.radian**2)
    system.addForce(angle)

    repulsive = openmm.CustomNonbondedForce("epsilon*(sigma/r)^12; sigma=1.0; epsilon=1.0")
    repulsive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    repulsive.setCutoffDistance(3.0*unit.nanometer)
    for _ in range(N):
        repulsive.addParticle([])
    system.addForce(repulsive)

    smc_force = openmm.HarmonicBondForce()
    max_bonds = n_lefs + 10
    for _ in range(max_bonds):
        smc_force.addBond(0, 1, 0.5*unit.nanometer, 0.0)
    system.addForce(smc_force)

    R = (3*N/(4*np.pi*POLYMER["density"]))**(1/3)
    confine = openmm.CustomExternalForce(
        f"k*max(0,r-R)^2; r=sqrt(x*x+y*y+z*z); R={R}; k=1.0")
    for i in range(N):
        confine.addParticle(i, [])
    system.addForce(confine)

    integrator = openmm.LangevinMiddleIntegrator(
        300*unit.kelvin, POLYMER["collision_rate"]/unit.picosecond, 0.01*unit.picosecond)

    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        props = {"CudaDeviceIndex": str(gpu)}
        context = openmm.Context(system, integrator, platform, props)
    except Exception:
        platform = openmm.Platform.getPlatformByName("CPU")
        context = openmm.Context(system, integrator, platform)

    context.setPositions(coords * unit.nanometer)
    context.setVelocitiesToTemperature(300*unit.kelvin)

    conformations = []
    smc_k = 100.0*unit.kilojoule_per_mole/unit.nanometer**2

    for block in range(n_blocks):
        lef_sim.step()
        bonds = lef_sim.get_bonds()

        for bi in range(max_bonds):
            if bi < len(bonds):
                l, r = bonds[bi]
                smc_force.setBondParameters(bi, l, r,
                    SMC_BOND["bond_dist"]*unit.nanometer, smc_k)
            else:
                smc_force.setBondParameters(bi, 0, 1, 0.5*unit.nanometer, 0.0)
        smc_force.updateParametersInContext(context)

        integrator.step(md_steps)

        if block % save_every == 0:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            conformations.append(pos.copy())

        if block % 500 == 0:
            logger.info(f"  block {block}/{n_blocks}, frames={len(conformations)}")

    with h5py.File(os.path.join(shard_dir, "conformations.h5"), "w") as hf:
        for i, conf in enumerate(conformations):
            hf.create_dataset(f"frame_{i}", data=conf)
        hf.attrs["N"] = N
        hf.attrs["n_frames"] = len(conformations)
    logger.info(f"  Saved {len(conformations)} frames")


def _run_shard_lef_only(lef_sim, N, n_blocks, save_every, shard_dir):
    """LEF-only fallback shard."""
    contact_map = np.zeros((N, N), dtype=np.float64)
    saved = 0

    for block in range(n_blocks):
        lef_sim.step()
        if block % save_every == 0:
            for (l, r) in lef_sim.get_bonds():
                contact_map[l, r] += 1
                contact_map[r, l] += 1
            saved += 1

        if block % 1000 == 0:
            logger.info(f"  block {block}/{n_blocks}, saved={saved}")

    np.save(os.path.join(shard_dir, "lef_contact_map.npy"), contact_map)
    logger.info(f"  Saved LEF contact map ({saved} frames)")


def main():
    parser = argparse.ArgumentParser(description="Sharded simulation runner")
    parser.add_argument("--condition", type=str, default=None,
                        help="Condition name from SIMULATION_CONDITIONS")
    parser.add_argument("--params", type=str, default=None,
                        help="Parameter set name (legacy, use --condition instead)")
    parser.add_argument("--ctcf-type", type=str, default=None,
                        choices=["mESC", "neuron"])
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--n-shards", type=int, required=True)
    parser.add_argument("--shard-index-offset", type=int, default=0,
                        help="Added to --shard-index for directory naming and seeding. "
                             "Use when running additional shards to extend an existing run "
                             "(e.g. previously ran --n-shards 2; now add 2 more with "
                             "--n-shards 2 --shard-index-offset 2).")
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--output", type=str, required=True)

    args = parser.parse_args()

    if not args.condition and not args.params:
        parser.error("Must provide either --condition or --params")
    params, ctcf_type = resolve_condition(
        condition_name=args.condition,
        params_name=args.params,
        ctcf_type=args.ctcf_type,
    )

    run_shard(params, args.replicate, args.shard_index, args.n_shards,
              args.gpu, args.output, ctcf_type=ctcf_type,
              shard_index_offset=args.shard_index_offset)


if __name__ == "__main__":
    main()
