#!/usr/bin/env python
"""
Main simulation script: run loop extrusion + 3D polymer simulation using polychrom.

This script:
  1. Initializes a polymer chain (representing the Sox2 locus, chr3:34-36 Mb)
  2. Runs loop extrusion dynamics (1D LEF model)
  3. Feeds LEF bonds into the 3D polymer simulation
  4. Saves conformations for downstream contact map analysis

Usage:
    python run_simulation.py --params mESC --replicate 0 --gpu 0
    python run_simulation.py --params CN_long_residency --replicate 0 --gpu 0

Requirements:
    pip install polychrom openmm h5py
"""

import os
import sys
import argparse
import time
import json
import logging
import numpy as np
import h5py

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import (
    POLYMER, SMC_BOND, SIM_RUN, N_MONOMERS, TILING,
    ALL_PARAM_SETS, SIMULATION_CONDITIONS,
    get_ctcf_arrays, get_tiled_ctcf_arrays, get_condition,
    CHROM, REGION_START, REGION_END, RESOLUTION,
)
from scripts.lef_dynamics import LEFSimulator

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def get_param_set(name: str) -> dict:
    """Look up a parameter set by name."""
    for ps in ALL_PARAM_SETS:
        if ps["name"] == name:
            return ps
    raise ValueError(f"Unknown parameter set: {name}. "
                     f"Available: {[p['name'] for p in ALL_PARAM_SETS]}")


def resolve_condition(condition_name: str = None, params_name: str = None,
                      ctcf_type: str = None):
    """
    Resolve simulation parameters and CTCF cell type.

    Can be called in two ways:
      1. --condition <name>  → looks up SIMULATION_CONDITIONS for both params + CTCF type
      2. --params <name> --ctcf-type <type>  → manual pairing (legacy / custom)

    Returns (params_dict, ctcf_type_str).
    """
    if condition_name:
        cond = get_condition(condition_name)
        return cond["params"], cond["ctcf_type"]
    elif params_name:
        params = get_param_set(params_name)
        ct = ctcf_type or "mESC"  # default to mESC for backward compat
        return params, ct
    else:
        raise ValueError("Must provide either --condition or --params")


def initialize_polymer(N: int, density: float = 0.2) -> np.ndarray:
    """
    Create an initial polymer configuration using a random walk
    confined to a sphere.

    Parameters
    ----------
    N : int
        Number of monomers.
    density : float
        Volume fraction.

    Returns
    -------
    coords : np.ndarray of shape (N, 3)
    """
    try:
        from polychrom.starting_conformations import grow_cubic
        return grow_cubic(N, int(N**0.5))
    except ImportError:
        logger.warning("polychrom not installed — using simple random walk initialization")
        # Simple confined random walk fallback
        R = (3 * N / (4 * np.pi * density)) ** (1.0 / 3)
        coords = np.zeros((N, 3))
        for i in range(1, N):
            step = np.random.randn(3)
            step /= np.linalg.norm(step)
            coords[i] = coords[i - 1] + step
            # Soft confinement
            dist = np.linalg.norm(coords[i])
            if dist > R:
                coords[i] *= R / dist
        return coords


def run_simulation_polychrom(params: dict, replicate: int, output_dir: str,
                             gpu: int = 0, ctcf_type: str = "mESC",
                             force_lef_only: bool = False,
                             resume: bool = False,
                             n_blocks_override: int | None = None):
    """
    Run full polychrom simulation with loop extrusion.

    Parameters
    ----------
    params : dict
        Cohesin parameter set (from parameters.py).
    replicate : int
        Replicate index.
    output_dir : str
        Directory to save results.
    gpu : int
        CUDA device index.
    ctcf_type : str
        Which CTCF site set to use: "mESC" or "neuron".
    force_lef_only : bool
        If True, skip the 3D polymer and only write the LEF contact map /
        trajectories. Used when ``--engine lef_only`` is requested so the
        polychrom path is never taken even if the library is installed.
    """
    if force_lef_only:
        HAS_POLYCHROM = False
        logger.info("force_lef_only=True — skipping 3D polymer, writing LEF data only")
    else:
        try:
            import polychrom
            from polychrom.simulation import Simulation
            from polychrom.forces import (
                harmonic_bonds,
                angle_force,
                polynomial_repulsive,
                spherical_confinement,
            )
            from polychrom.hdf5_format import HDF5Reporter
            HAS_POLYCHROM = True
        except ImportError:
            HAS_POLYCHROM = False
            logger.warning("polychrom not available — running LEF-only simulation")

    # --- Tiling: simulate a large chromosome with N_TILES copies of the locus ---
    # Numeric values below are the current defaults from configs/parameters.py
    # (TILING dict); they are kept as comments only to make the intent clear.
    # The truth is always TILING at runtime.
    N_locus = POLYMER["N"]                      # 2000 (original locus)
    N = TILING["chrom_size"]                     # 70000 (tiled chromosome)
    n_tiles = TILING["n_tiles"]                  # 28
    tile_size = TILING["tile_size"]              # 2500
    pad = TILING["padding"]                      # 250

    name = params["name"]
    # Include CTCF type in the directory name so mESC-CTCF and neuron-CTCF
    # runs are stored separately even when using the same cohesin parameters
    run_dir = os.path.join(output_dir, f"{name}_ctcf-{ctcf_type}_rep{replicate}")
    os.makedirs(run_dir, exist_ok=True)

    # Save parameters (including which CTCF set was used and tiling info)
    with open(os.path.join(run_dir, "params.json"), "w") as f:
        json.dump({**params, "ctcf_type": ctcf_type,
                    "polymer": POLYMER, "sim_run": SIM_RUN, "tiling": TILING,
                    "replicate": replicate}, f, indent=2, default=str)

    # --- Resume detection ---
    # Resume kicks in only if (a) the caller asked for it AND (b) there is
    # already at least one blocks_*.h5 in run_dir to resume from. Without
    # both, we fall through to a fresh run with HDF5Reporter overwrite=True.
    import glob as _glob
    existing_block_files = sorted(_glob.glob(os.path.join(run_dir, "blocks_*.h5")))
    should_resume = bool(resume and existing_block_files)
    if resume and not should_resume:
        logger.info(f"--resume requested, but {run_dir} has no blocks_*.h5 yet; fresh start.")
    if should_resume:
        logger.info(f"--resume: {len(existing_block_files)} block file(s) already present in {run_dir}")

    # --- Per-rep deterministic seeding ---
    # Seed every stochastic component with a function of `replicate` so that
    # (a) different reps are guaranteed statistically independent, and
    # (b) any individual rep is reproducible if re-run with the same number.
    # The LEFSimulator already takes an explicit rng_seed; we additionally
    # seed numpy (used by polychrom.starting_conformations.grow_cubic and any
    # other np.random.* call) and the OpenMM integrator (Langevin noise) below
    # right after Simulation() construction.
    rep_seed = 42 + replicate * 1000
    np.random.seed(rep_seed)

    # --- Initialize LEF simulator on the FULL tiled chromosome ---
    ctcf_pos, ctcf_ori = get_tiled_ctcf_arrays(cell_type=ctcf_type)
    n_lefs = max(1, N // params["separation"])

    lef_sim = LEFSimulator(
        N=N,
        n_lefs=n_lefs,
        lifetime=params["lifetime"],
        ctcf_positions=ctcf_pos,
        ctcf_orientations=ctcf_ori,
        ctcf_capture=params["ctcf_capture"],
        ctcf_release=params["ctcf_release"],
        rng_seed=rep_seed,
    )

    logger.info(f"LEF simulator initialized: {n_lefs} LEFs on {N}-monomer tiled chromosome "
                f"({n_tiles} tiles × {tile_size} monomers), "
                f"lifetime={params['lifetime']}, separation={params['separation']}, "
                f"CTCF={ctcf_type} ({len(ctcf_pos)} sites total, "
                f"{len(ctcf_pos)//n_tiles} per tile)")

    # --- Warm up LEFs ---
    logger.info("Warming up LEF dynamics...")
    for _ in range(SIM_RUN["warmup_blocks"]):
        for __ in range(SIM_RUN["smc_steps_per_block"]):
            lef_sim.step()

    if HAS_POLYCHROM:
        # =====================================================================
        # FULL POLYCHROM SIMULATION
        # Pattern adapted from Hansen lab: Yang JH, Brandão HB, Hansen AS
        # (2023) Nat Commun 14:1913 "DNA double-strand break end synapsis by
        # DNA loop extrusion", DOI 10.1038/s41467-023-37583-w. See
        # 3D_PolymerSimulationCode(withLoopExtrusion).py in
        # github.com/ahansenlab/DNA_break_synapsis_models.
        # =====================================================================
        logger.info("Starting polychrom 3D simulation...")

        total_blocks = n_blocks_override if n_blocks_override is not None else SIM_RUN["total_blocks"]
        if n_blocks_override is not None:
            logger.info(f"total_blocks overridden via --n-blocks: {total_blocks} "
                        f"(SIM_RUN default would be {SIM_RUN['total_blocks']})")
        save_every = SIM_RUN["save_every"]
        smc_steps = SIM_RUN["smc_steps_per_block"]
        md_steps = SIM_RUN["md_steps_per_block"]

        # --- Pre-compute all LEF bond states ---
        # LEF dynamics are 1D and independent of 3D positions, so we can
        # pre-compute all bond configurations before starting the 3D sim.
        logger.info("Pre-computing LEF bond states for all blocks...")
        all_bond_sets = []
        for _ in range(total_blocks):
            for __ in range(smc_steps):
                lef_sim.step()
            bonds = lef_sim.get_bonds()
            all_bond_sets.append([(int(l), int(r)) for l, r in bonds])

        # Find all unique bonds across the entire simulation
        unique_bonds = list(set(sum(all_bond_sets, [])))
        logger.info(f"Pre-computed {total_blocks} bond states, "
                     f"{len(unique_bonds)} unique bond pairs")

        # --- Initialize polymer and simulation ---
        # On resume, load the last saved conformation from blocks_*.h5 instead
        # of growing a fresh polymer; the HDF5Reporter is opened in append
        # mode (overwrite=False, check_exists=False) and `continue_trajectory`
        # gives us back (last_block_idx, dict_with_pos) so the next saved
        # block lands at last_block_idx+1.
        if should_resume:
            reporter = HDF5Reporter(
                folder=run_dir, max_data_length=100,
                overwrite=False, check_exists=False, blocks_only=True,
            )
            last_block_idx, last_data = reporter.continue_trajectory()
            start_block = int(last_block_idx) + 1
            coords = np.asarray(last_data["pos"])
            logger.info(f"Resumed at block {start_block} "
                        f"(last saved block was {last_block_idx}, "
                        f"loaded conformation shape {coords.shape})")
        else:
            reporter = HDF5Reporter(
                folder=run_dir, max_data_length=100,
                overwrite=True, blocks_only=True,
            )
            start_block = 0
            coords = initialize_polymer(N, POLYMER["density"])

        sim = Simulation(
            platform=POLYMER["platform"],
            GPU=str(gpu),
            N=N,
            error_tol=POLYMER["error_tol"],
            collision_rate=POLYMER["collision_rate"],
            integrator="variableLangevin",
            reporters=[reporter],
            PBCbox=False,
            precision="mixed",
        )

        # Seed the OpenMM integrator's Langevin noise deterministically per
        # replicate so the thermal trajectory is also reproducible.
        try:
            sim.integrator.setRandomNumberSeed(rep_seed)
            logger.info(f"OpenMM integrator seed set to {rep_seed}")
        except Exception as e:
            logger.warning(f"Could not set integrator seed: {e}")

        sim.set_data(coords)

        # --- Add forces ---
        sim.add_force(spherical_confinement(
            sim, density=POLYMER["density"], k=1.0,
        ))

        from polychrom import forcekits
        sim.add_force(
            forcekits.polymer_chains(
                sim,
                chains=[(0, N, False)],  # one linear chain (tiled chromosome)
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
        # Use polychrom's own scaling factors for correct unit handling
        kbond = sim.kbondScalingFactor / (SMC_BOND["wiggle_dist"] ** 2)
        bondDist = SMC_BOND["bond_dist"] * sim.length_scale
        activeParams = {"length": bondDist, "k": kbond}
        inactiveParams = {"length": bondDist, "k": 0}

        bondForce = sim.force_dict["harmonic_bonds"]
        # On resume, initialise bonds at the state they had during the LAST
        # already-saved block, so the diff update inside the main loop only
        # carries the new transitions from start_block onward.
        bonds_at_start = all_bond_sets[start_block - 1] if should_resume else all_bond_sets[0]
        bond_to_ind = {}
        for bond in unique_bonds:
            paramset = activeParams if (bond in bonds_at_start) else inactiveParams
            ind = bondForce.addBond(bond[0], bond[1], **paramset)
            bond_to_ind[bond] = ind

        # --- Create the OpenMM context ---
        # Skip energy minimisation on resume — the loaded conformation is
        # already a relaxed post-MD state.
        if not should_resume:
            sim.local_energy_minimization()
            logger.info("Energy minimization complete, starting production run...")
        else:
            logger.info(f"Skipping energy minimization (resume from block {start_block}).")

        # --- Main simulation loop ---
        saved_count = 0
        cur_bonds = bonds_at_start
        # Carry the bond set across loop iterations so we only pay one set()
        # conversion per block (the bonds-add/remove diff is set arithmetic
        # below; see the per-block comment for the rationale).
        cur_bonds_set = set(cur_bonds)

        # Optional per-loop-section timing breakdown. Enable via env var
        # PROFILE_BLOCK_LOOP=1; off by default (production runs untouched).
        # Set PROFILE_REPORT_EVERY (default 50) to control the report cadence.
        import time as _time
        profile_loop = os.environ.get("PROFILE_BLOCK_LOOP", "0") == "1"
        report_every = int(os.environ.get("PROFILE_REPORT_EVERY", "50"))
        t_diff = t_set = t_update = t_md = 0.0
        n_add_total = n_rem_total = 0

        for block in range(start_block, total_blocks):
            new_bonds = all_bond_sets[block]

            if profile_loop:
                t0 = _time.perf_counter()

            # Toggle changed bonds only (efficient delta update).
            # Set-based diff is O(n) vs the previous O(n²) list-comprehension
            # membership tests. cur_bonds_set is maintained across iterations
            # so we don't pay set() conversion twice per loop.
            new_bonds_set = set(new_bonds)
            bonds_remove = cur_bonds_set - new_bonds_set
            bonds_add = new_bonds_set - cur_bonds_set
            cur_bonds_set = new_bonds_set

            if profile_loop:
                t1 = _time.perf_counter(); t_diff += t1 - t0
                n_add_total += len(bonds_add); n_rem_total += len(bonds_remove)

            for bond in bonds_add:
                ind = bond_to_ind[bond]
                bondForce.setBondParameters(ind, bond[0], bond[1], **activeParams)
            for bond in bonds_remove:
                ind = bond_to_ind[bond]
                bondForce.setBondParameters(ind, bond[0], bond[1], **inactiveParams)

            if profile_loop:
                t2 = _time.perf_counter(); t_set += t2 - t1

            bondForce.updateParametersInContext(sim.context)

            if profile_loop:
                t3 = _time.perf_counter(); t_update += t3 - t2

            # cur_bonds_set is already updated above (= new_bonds_set);
            # cur_bonds list kept for any non-loop consumer that may read it.
            cur_bonds = new_bonds

            # Run MD and save
            sim.do_block(md_steps)
            saved_count += 1

            if profile_loop:
                t4 = _time.perf_counter(); t_md += t4 - t3

            if profile_loop and (block + 1) % report_every == 0:
                n = report_every
                total = t_diff + t_set + t_update + t_md
                logger.info(
                    f"[profile] last {n} blocks: "
                    f"diff={1000*t_diff/n:.1f} ms, "
                    f"setBondParameters={1000*t_set/n:.1f} ms "
                    f"({n_add_total + n_rem_total} bond-toggles), "
                    f"updateParametersInContext={1000*t_update/n:.1f} ms, "
                    f"do_block(MD,{md_steps})={1000*t_md/n:.1f} ms | "
                    f"per-block total={1000*total/n:.1f} ms; "
                    f"shares: diff {100*t_diff/total:.0f}% / "
                    f"setBP {100*t_set/total:.0f}% / "
                    f"updateP {100*t_update/total:.0f}% / "
                    f"MD {100*t_md/total:.0f}%"
                )
                t_diff = t_set = t_update = t_md = 0.0
                n_add_total = n_rem_total = 0

            if block % 500 == 0:
                logger.info(f"Block {block}/{total_blocks} | "
                            f"LEFs: {len(new_bonds)} bonds | saved: {saved_count}")

        reporter.dump_data()
        logger.info(f"Polychrom simulation complete. Saved {saved_count} conformations.")

    else:
        # =====================================================================
        # LEF-ONLY SIMULATION (no polychrom — generate contact maps from 1D)
        # =====================================================================
        logger.info("Running LEF-only simulation (no 3D polymer)...")

        total_blocks = n_blocks_override if n_blocks_override is not None else SIM_RUN["total_blocks"]
        save_every = SIM_RUN["save_every"]
        warmup = SIM_RUN["warmup_blocks"]

        # We'll accumulate a contact-like matrix from LEF bridging
        contact_map = np.zeros((N, N), dtype=np.float64)
        saved_count = 0
        all_conformations = []

        for block in range(total_blocks):
            for _ in range(SIM_RUN["smc_steps_per_block"]):
                lef_sim.step()

            if block >= warmup and block % save_every == 0:
                bonds = lef_sim.get_bonds()
                for (l, r) in bonds:
                    contact_map[l, r] += 1
                    contact_map[r, l] += 1

                # Also store loop sizes for analysis
                all_conformations.append({
                    "block": block,
                    "bonds": bonds,
                    "loop_sizes": lef_sim.get_loop_sizes().tolist(),
                })
                saved_count += 1

            if block % 1000 == 0:
                logger.info(f"Block {block}/{total_blocks} | saved: {saved_count}")

        # Save results
        np.save(os.path.join(run_dir, "lef_contact_map.npy"), contact_map)

        with h5py.File(os.path.join(run_dir, "lef_trajectories.h5"), "w") as hf:
            for i, conf in enumerate(all_conformations):
                grp = hf.create_group(f"frame_{i}")
                bonds_arr = np.array(conf["bonds"]) if conf["bonds"] else np.array([]).reshape(0, 2)
                grp.create_dataset("bonds", data=bonds_arr)
                grp.create_dataset("loop_sizes", data=np.array(conf["loop_sizes"]))
                grp.attrs["block"] = conf["block"]

        logger.info(f"LEF-only simulation complete. Saved {saved_count} frames.")

    logger.info(f"Results saved to: {run_dir}")
    return run_dir


def run_simulation_openmm_standalone(params: dict, replicate: int, output_dir: str,
                                     gpu: int = 0, ctcf_type: str = "mESC",
                                     resume: bool = False,
                                     n_blocks_override: int | None = None):
    """
    Standalone OpenMM simulation without polychrom dependency.
    Uses raw OpenMM API for the 3D polymer, with LEF bonds updated dynamically.

    This is the recommended fallback if polychrom is not installed.
    Note: this path writes a single conformations.h5 file (not per-block),
    so resume is not implemented here — the flag is accepted but ignored
    with a warning, so callers don't need to branch.
    """
    if resume:
        logger.warning("resume=True ignored: standalone OpenMM path writes "
                       "a single conformations.h5 and has no per-block resume.")
    try:
        import openmm
        import openmm.app as app
        import openmm.unit as unit
        HAS_OPENMM = True
    except ImportError:
        HAS_OPENMM = False

    if not HAS_OPENMM:
        logger.info("OpenMM not available. Falling back to LEF-only simulation.")
        return run_simulation_polychrom(params, replicate, output_dir, gpu,
                                        ctcf_type=ctcf_type, resume=resume,
                                        n_blocks_override=n_blocks_override)

    N = POLYMER["N"]
    name = params["name"]
    run_dir = os.path.join(output_dir, f"{name}_ctcf-{ctcf_type}_rep{replicate}")
    os.makedirs(run_dir, exist_ok=True)

    with open(os.path.join(run_dir, "params.json"), "w") as f:
        json.dump({**params, "ctcf_type": ctcf_type,
                    "polymer": POLYMER, "sim_run": SIM_RUN,
                    "replicate": replicate}, f, indent=2, default=str)

    # Per-rep deterministic seeding (mirrors run_simulation_polychrom).
    rep_seed = 42 + replicate * 1000
    np.random.seed(rep_seed)

    # Initialize LEFs
    ctcf_pos, ctcf_ori = get_ctcf_arrays(cell_type=ctcf_type)
    n_lefs = max(1, N // params["separation"])
    lef_sim = LEFSimulator(
        N=N, n_lefs=n_lefs, lifetime=params["lifetime"],
        ctcf_positions=ctcf_pos, ctcf_orientations=ctcf_ori,
        ctcf_capture=params["ctcf_capture"], ctcf_release=params["ctcf_release"],
        rng_seed=rep_seed,
    )

    # Warm up LEFs
    for _ in range(SIM_RUN["warmup_blocks"]):
        lef_sim.step()

    # --- Build OpenMM system ---
    system = openmm.System()

    # Add particles
    mass = 100.0 * unit.amu
    for _ in range(N):
        system.addParticle(mass)

    # Backbone harmonic bonds
    backbone = openmm.HarmonicBondForce()
    backbone.setUsesPeriodicBoundaryConditions(False)
    for i in range(N - 1):
        backbone.addBond(i, i + 1,
                         POLYMER["bond_length"] * unit.nanometer,
                         300.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(backbone)

    # Angle force for stiffness
    angle = openmm.HarmonicAngleForce()
    for i in range(N - 2):
        angle.addAngle(i, i + 1, i + 2,
                       np.pi * unit.radian,
                       POLYMER["angle_force"] * unit.kilojoule_per_mole / unit.radian**2)
    system.addForce(angle)

    # Excluded volume (soft repulsion via custom nonbonded)
    repulsive = openmm.CustomNonbondedForce(
        "epsilon * (sigma/r)^12; sigma=1.0; epsilon=1.0"
    )
    repulsive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    repulsive.setCutoffDistance(3.0 * unit.nanometer)
    for _ in range(N):
        repulsive.addParticle([])
    system.addForce(repulsive)

    # SMC bond force (will be updated dynamically)
    smc_force = openmm.HarmonicBondForce()
    smc_force.setUsesPeriodicBoundaryConditions(False)
    # Pre-allocate bond slots
    max_bonds = n_lefs + 10
    for _ in range(max_bonds):
        smc_force.addBond(0, 1, 0.5 * unit.nanometer, 0.0)  # inactive
    system.addForce(smc_force)

    # Spherical confinement
    R_confine = (3 * N / (4 * np.pi * POLYMER["density"])) ** (1.0 / 3)
    confine = openmm.CustomExternalForce(
        f"k * max(0, r - R)^2; r = sqrt(x*x + y*y + z*z); R = {R_confine}; k = 1.0"
    )
    for i in range(N):
        confine.addParticle(i, [])
    system.addForce(confine)

    # Integrator
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        POLYMER["collision_rate"] / unit.picosecond,
        0.01 * unit.picosecond,
    )
    integrator.setRandomNumberSeed(rep_seed)

    # Platform
    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        properties = {"CudaDeviceIndex": str(gpu)}
        context = openmm.Context(system, integrator, platform, properties)
    except Exception:
        platform = openmm.Platform.getPlatformByName("CPU")
        context = openmm.Context(system, integrator, platform)

    # Set initial positions
    coords = initialize_polymer(N, POLYMER["density"])
    context.setPositions(coords * unit.nanometer)
    context.setVelocitiesToTemperature(300 * unit.kelvin)

    # --- Main simulation loop ---
    total_blocks = n_blocks_override if n_blocks_override is not None else SIM_RUN["total_blocks"]
    save_every = SIM_RUN["save_every"]
    warmup = SIM_RUN["warmup_blocks"]
    md_steps = SIM_RUN["md_steps_per_block"]

    conformations = []
    saved_count = 0

    logger.info(f"Starting OpenMM simulation: {name}, replicate {replicate}")

    for block in range(total_blocks):
        # Update LEF bonds
        lef_sim.step()
        bonds = lef_sim.get_bonds()

        # Update SMC force bonds
        smc_k = 100.0 * unit.kilojoule_per_mole / unit.nanometer**2
        for bi in range(max_bonds):
            if bi < len(bonds):
                l, r = bonds[bi]
                smc_force.setBondParameters(bi, l, r,
                                            SMC_BOND["bond_dist"] * unit.nanometer,
                                            smc_k)
            else:
                smc_force.setBondParameters(bi, 0, 1,
                                            0.5 * unit.nanometer,
                                            0.0)  # zero force
        smc_force.updateParametersInContext(context)

        # Run MD
        integrator.step(md_steps)

        # Save conformations
        if block >= warmup and block % save_every == 0:
            state = context.getState(getPositions=True)
            pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            conformations.append(pos.copy())
            saved_count += 1

        if block % 500 == 0:
            logger.info(f"Block {block}/{total_blocks} | bonds: {len(bonds)} | saved: {saved_count}")

    # Save conformations
    with h5py.File(os.path.join(run_dir, "conformations.h5"), "w") as hf:
        for i, conf in enumerate(conformations):
            hf.create_dataset(f"frame_{i}", data=conf)
        hf.attrs["N"] = N
        hf.attrs["n_frames"] = len(conformations)
        hf.attrs["params"] = json.dumps(params)

    logger.info(f"Simulation complete. Saved {saved_count} conformations to {run_dir}")
    return run_dir


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Run cohesin loop extrusion simulation",
        epilog="Examples:\n"
               "  # Using a predefined condition (recommended):\n"
               "  python run_simulation.py --condition mESC_ctrl --replicate 0 --gpu 0\n"
               "  python run_simulation.py --condition CN_long_residency_neuron_ctcf --replicate 0\n\n"
               "  # Legacy: manual pairing of params + CTCF type:\n"
               "  python run_simulation.py --params mESC --ctcf-type neuron --replicate 0\n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # New condition-based interface (preferred)
    parser.add_argument("--condition", type=str, default=None,
                        help="Simulation condition name from SIMULATION_CONDITIONS "
                             "(e.g., mESC_ctrl, CN_long_residency_neuron_ctcf). "
                             "This sets both the cohesin parameters and the CTCF site set.")

    # Legacy params-based interface (still supported)
    parser.add_argument("--params", type=str, default=None,
                        help="Parameter set name (e.g., mESC, CN_long_residency). "
                             "Use with --ctcf-type for manual pairing.")
    parser.add_argument("--ctcf-type", type=str, default=None,
                        choices=["mESC", "neuron"],
                        help="Which CTCF site set to use (default: mESC). "
                             "Only used with --params, ignored with --condition.")

    parser.add_argument("--replicate", type=int, default=0,
                        help="Replicate index")
    parser.add_argument("--gpu", type=int, default=0,
                        help="GPU device index")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory (default: ../results/)")
    parser.add_argument("--resume", action="store_true",
                        help="If the rep dir already contains blocks_*.h5, "
                             "continue from the last saved block instead of "
                             "wiping. No-op when the rep dir is fresh. "
                             "Required for the catch-up walltime-restart loop.")
    parser.add_argument("--n-blocks", type=int, default=None,
                        help="Override SIM_RUN['total_blocks'] for this run. "
                             "Used by the catch-up wrapper to size each rep "
                             "to fit comfortably in the 48h walltime "
                             "(default 700 blocks ≈ 22h on V100, "
                             "yielding 70 000 saved conformations).")
    parser.add_argument("--engine", type=str, default="auto",
                        choices=["auto", "polychrom", "openmm", "lef_only"],
                        help="Simulation engine to use")

    args = parser.parse_args()

    # Resolve condition → (params, ctcf_type)
    if not args.condition and not args.params:
        parser.error("Must provide either --condition or --params")
    params, ctcf_type = resolve_condition(
        condition_name=args.condition,
        params_name=args.params,
        ctcf_type=args.ctcf_type,
    )

    if args.output is None:
        args.output = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    os.makedirs(args.output, exist_ok=True)

    logger.info(f"Running simulation: {params['name']} + CTCF={ctcf_type} (rep {args.replicate})")
    logger.info(f"  Lifetime: {params['lifetime']}, Separation: {params['separation']}")
    logger.info(f"  CTCF capture: {params['ctcf_capture']}, release: {params['ctcf_release']}")

    run_kwargs = dict(params=params, replicate=args.replicate,
                      output_dir=args.output, gpu=args.gpu, ctcf_type=ctcf_type,
                      resume=args.resume, n_blocks_override=args.n_blocks)

    if args.engine == "lef_only":
        # Skip the 3D polymer even if polychrom is installed: produce only
        # the 1D LEF contact map and bond trajectories.
        run_simulation_polychrom(**run_kwargs, force_lef_only=True)
    elif args.engine == "polychrom":
        run_simulation_polychrom(**run_kwargs)
    elif args.engine == "openmm":
        run_simulation_openmm_standalone(**run_kwargs)
    else:
        # Auto-detect best engine
        try:
            import polychrom
            run_simulation_polychrom(**run_kwargs)
        except ImportError:
            try:
                import openmm
                run_simulation_openmm_standalone(**run_kwargs)
            except ImportError:
                logger.warning("Neither polychrom nor openmm available. Running LEF-only.")
                run_simulation_polychrom(**run_kwargs, force_lef_only=True)


if __name__ == "__main__":
    main()
