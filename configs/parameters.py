"""
Simulation parameters loader — reads config/SimConfig.yaml.

Presents the same public API as the original hardcoded parameters.py so all
downstream scripts work without modification:

    ACTIVE_LOCUS, LOCI, CHROM, REGION_START, REGION_END, RESOLUTION, N_MONOMERS
    POLYMER, SMC_BOND, TILING, SIM_RUN
    SIMULATION_CONDITIONS, ALL_PARAM_SETS
    CTCF_SITES_CONTROL, CTCF_SITES_SORBITOL, CTCF_SITES
    CTCF_BED_CONTROL, CTCF_BED_SORBITOL, CTCF_BED_PROMOTER
    get_condition(), list_conditions()
    get_ctcf_arrays(), load_ctcf_from_bed(), load_promoter_sites()
    get_tiled_ctcf_arrays()
    DATA, MSD_PROBE, MSD_FIT, CALIBRATION
"""

import numpy as np
import os as _os

# ── Locate and load SimConfig.yaml ────────────────────────────────────────────
_REPO_ROOT   = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), ".."))
_CONFIG_PATH = _os.path.join(_REPO_ROOT, "config", "SimConfig.yaml")

try:
    import yaml as _yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml (or conda install pyyaml)")

with open(_CONFIG_PATH) as _f:
    _cfg = _yaml.safe_load(_f)


# =============================================================================
# LOCUS
# =============================================================================
_loc         = _cfg["locus"]
ACTIVE_LOCUS = _loc["name"]

LOCI = {
    _loc["name"]: dict(
        chrom        = _loc["chrom"],
        region_start = int(_loc["region_start"]),
        region_end   = int(_loc["region_end"]),
        label        = _loc.get("label", ""),
        anchor_bp    = int(_loc.get("anchor_bp",
                           (_loc["region_start"] + _loc["region_end"]) // 2)),
    )
}

_locus       = LOCI[ACTIVE_LOCUS]
CHROM        = _locus["chrom"]
REGION_START = _locus["region_start"]
REGION_END   = _locus["region_end"]
RESOLUTION   = int(_cfg.get("resolution", 1000))
N_MONOMERS   = (REGION_END - REGION_START) // RESOLUTION

_ASSEMBLY    = _cfg.get("genome_assembly", "hg38")


# =============================================================================
# POLYMER PHYSICS
# =============================================================================
_poly = _cfg.get("polymer", {})
POLYMER = dict(
    N                = N_MONOMERS,
    density          = _poly.get("density",          0.2),
    bond_length      = _poly.get("bond_length",      1.0),
    bond_wiggle_dist = _poly.get("bond_wiggle_dist", 0.05),
    angle_force      = _poly.get("angle_force",      1.5),
    collision_rate   = _poly.get("collision_rate",   0.03),
    error_tol        = _poly.get("error_tol",        0.003),
    platform         = _poly.get("platform",         "CUDA"),
)

_smc = _cfg.get("smc_bond", {})
SMC_BOND = dict(
    wiggle_dist = _smc.get("wiggle_dist", 0.2),
    bond_dist   = _smc.get("bond_dist",   0.5),
)


# =============================================================================
# TILING
# =============================================================================
_til = _cfg.get("tiling", {})
TILING = dict(
    tile_size  = int(_til.get("tile_size",  2500)),
    chrom_size = int(_til.get("chrom_size", 70000)),
    padding    = int(_til.get("padding",    250)),
)
TILING["n_tiles"] = TILING["chrom_size"] // TILING["tile_size"]

assert TILING["tile_size"] == N_MONOMERS + 2 * TILING["padding"], (
    f"tile_size ({TILING['tile_size']}) must equal "
    f"N_MONOMERS ({N_MONOMERS}) + 2*padding ({2*TILING['padding']}). "
    f"Check locus coordinates and tiling section in SimConfig.yaml."
)


# =============================================================================
# SIMULATION RUN PARAMETERS
# =============================================================================
_run = _cfg.get("sim_run", {})
SIM_RUN = dict(
    smc_steps_per_block  = int(_run.get("smc_steps_per_block",  1)),
    md_steps_per_block   = int(_run.get("md_steps_per_block",   250)),
    total_blocks         = int(_run.get("total_blocks",         50000)),
    warmup_blocks        = int(_run.get("warmup_blocks",        4000)),
    save_every           = int(_run.get("save_every",           10)),
    restart_milker_every = int(_run.get("restart_milker_every", 100)),
    contact_radius       = float(_run.get("contact_radius",     3.0)),
    n_replicates         = int(_run.get("n_replicates", _cfg.get("n_replicates", 3))),
)


# =============================================================================
# CTCF BED FILE PATHS
# =============================================================================
# Lives in data/ctcf_beds/; naming: ctcf_oriented_{assembly}_{ctcf_type}_{chrom}_{start}_{end}.bed
def _ctcf_bed_path(ctcf_type):
    return (
        f"data/ctcf_beds/ctcf_oriented_{_ASSEMBLY}_{ctcf_type}_"
        f"{CHROM}_{REGION_START}_{REGION_END}.bed"
    )

CTCF_BED_CONTROL  = _ctcf_bed_path("control")
CTCF_BED_SORBITOL = _ctcf_bed_path("sorbitol")
CTCF_BED_PROMOTER = (
    f"data/ctcf_beds/promoter_stall_sites_{CHROM}_{REGION_START}_{REGION_END}.bed"
)
AUTO_LOAD_CTCF_FROM_BED = True


# =============================================================================
# CTCF SITE PLACEHOLDERS (replaced at import time by autoload below)
# =============================================================================
CTCF_SITES_PLACEHOLDER = True
CTCF_SITES_CONTROL     = [(500, +1), (1000, -1), (1500, +1)]
CTCF_SITES_SORBITOL    = [(500, +1), (1500, +1)]
CTCF_SITES             = CTCF_SITES_CONTROL
PROMOTER_STALL_SITES   = []


# =============================================================================
# SIMULATION CONDITIONS (built from YAML)
# =============================================================================
SIMULATION_CONDITIONS = []

for _c in _cfg.get("conditions", []):
    _params = dict(
        name         = _c.get("params_name", _c["name"]),
        description  = _c.get("description", _c["name"]),
        lifetime     = int(_c["lifetime"]),
        separation   = int(_c["separation"]),
        ctcf_capture = float(_c["ctcf_capture"]),
        ctcf_release = float(_c["ctcf_release"]),
    )
    _cond = dict(
        name        = _c["name"],
        description = _c.get("description", _c["name"]),
        params      = _params,
        ctcf_type   = _c.get("ctcf_type", "control"),
        compare_hic = _c.get("compare_hic", None),
    )
    if "ctcf_extra" in _c:
        _cond["ctcf_extra"] = _c["ctcf_extra"]
    SIMULATION_CONDITIONS.append(_cond)

ALL_PARAM_SETS = [c["params"] for c in SIMULATION_CONDITIONS]


def get_condition(name):
    """Look up a simulation condition by name."""
    for cond in SIMULATION_CONDITIONS:
        if cond["name"] == name:
            return cond
    available = [c["name"] for c in SIMULATION_CONDITIONS]
    raise ValueError(f"Unknown condition '{name}'. Available: {available}")


def list_conditions():
    """Print all simulation conditions as a formatted table."""
    print(f"{'Condition':<45} {'CTCF':>10} {'Hi-C':>8}   Description")
    print("-" * 110)
    for c in SIMULATION_CONDITIONS:
        hic   = c["compare_hic"] or "—"
        extra = f"+{c['ctcf_extra']}" if c.get("ctcf_extra") else ""
        print(
            f"{c['name']:<45} "
            f"{c['ctcf_type'] + extra:>10} "
            f"{hic:>8}   "
            f"{c['description']}"
        )


# =============================================================================
# CTCF LOADING FUNCTIONS
# =============================================================================

def get_ctcf_arrays(condition="control"):
    """Return (positions, orientations) arrays for the given condition."""
    if CTCF_SITES_PLACEHOLDER:
        import warnings
        warnings.warn(
            f"Using PLACEHOLDER CTCF sites for condition='{condition}'! "
            "Run setup_data.sh to generate oriented BED files.",
            UserWarning,
        )
    sites = CTCF_SITES_CONTROL if condition == "control" else CTCF_SITES_SORBITOL
    positions    = np.array([s[0] for s in sites], dtype=int)
    orientations = np.array([s[1] for s in sites], dtype=int)
    return positions, orientations


def load_ctcf_from_bed(bed_path, condition="control",
                       region_start=None, resolution=None):
    """
    Load oriented CTCF sites from a BED6 file.
    Populates CTCF_SITES_CONTROL or CTCF_SITES_SORBITOL (and CTCF_SITES alias).
    """
    global CTCF_SITES_CONTROL, CTCF_SITES_SORBITOL, CTCF_SITES, CTCF_SITES_PLACEHOLDER

    if region_start is None:
        region_start = REGION_START
    if resolution is None:
        resolution = RESOLUTION

    sites = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 3 or fields[0] != CHROM:
                continue
            start  = int(fields[1])
            end    = int(fields[2])
            strand = fields[5] if len(fields) > 5 else "+"
            summit  = (start + end) // 2
            monomer = (summit - region_start) // resolution
            if 0 <= monomer < N_MONOMERS:
                sites.append((monomer, +1 if strand == "+" else -1))

    sites.sort(key=lambda x: x[0])

    if condition == "sorbitol":
        CTCF_SITES_SORBITOL = sites
    else:
        CTCF_SITES_CONTROL = sites
        CTCF_SITES = sites

    CTCF_SITES_PLACEHOLDER = False
    return sites


def load_promoter_sites(bed_path, region_start=None, resolution=None):
    """Load promoter-based cohesin stall sites (non-directional, orientation +1)."""
    global PROMOTER_STALL_SITES
    if region_start is None:
        region_start = REGION_START
    if resolution is None:
        resolution = RESOLUTION
    sites = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3 or fields[0] != CHROM:
                continue
            summit  = (int(fields[1]) + int(fields[2])) // 2
            monomer = (summit - region_start) // resolution
            strand  = fields[5] if len(fields) > 5 else "+"
            if 0 <= monomer < N_MONOMERS:
                sites.append((monomer, +1 if strand == "+" else -1))
    sites.sort(key=lambda x: x[0])
    PROMOTER_STALL_SITES = sites
    return sites


def get_tiled_ctcf_arrays(condition="control"):
    """Return CTCF positions/orientations tiled across the full simulated chromosome."""
    positions, orientations = get_ctcf_arrays(condition)
    pad       = TILING["padding"]
    tile_size = TILING["tile_size"]
    n_tiles   = TILING["n_tiles"]
    tiled_pos, tiled_ori = [], []
    for tile in range(n_tiles):
        offset = tile * tile_size + pad
        tiled_pos.extend((positions + offset).tolist())
        tiled_ori.extend(orientations.tolist())
    return np.array(tiled_pos, dtype=int), np.array(tiled_ori, dtype=int)


# =============================================================================
# HI-C DATA PATHS
# =============================================================================
DATA = dict(
    hic_control    = _cfg.get("hic_control",   ""),
    hic_sorbitol   = _cfg.get("hic_sorbitol",  ""),
    mcool_control  = "data/mcool/control.mcool",
    mcool_sorbitol = "data/mcool/sorbitol.mcool",
    resolution     = RESOLUTION,
)


# =============================================================================
# MSD + CALIBRATION (analysis-side; kept for backwards compatibility)
# =============================================================================
MSD_PROBE = dict(
    pairs_control   = None,
    pairs_sorbitol  = None,
    labels_control  = None,
    labels_sorbitol = None,
    auto = dict(
        anchor_bp              = _locus["anchor_bp"],
        min_sep_monomers       = 200,
        max_sep_monomers       = 1500,
        require_convergent     = True,
        n_pairs                = 2,
        require_conserved      = True,
        conserved_tol_monomers = 0,
        prefer_flanking        = True,
        disjoint_anchors       = True,
        relax_if_missing       = True,
    ),
)

MSD_FIT = dict(
    lag_min_frames     = 1,
    lag_max_frac       = 1 / 3,
    fit_lag_min        = 5,
    fit_lag_max_frac   = 0.25,
    min_n_lags_for_fit = 6,
)

CALIBRATION = dict(
    default_anchor = "hic",
    hic = dict(
        s_ref_bp             = 100_000,
        persistence_nm       = 50.0,
        nm_per_bp            = 0.34,
        fixed_bp_per_monomer = RESOLUTION,
    ),
    msd = dict(
        expt_msd_control  = None,
        expt_msd_sorbitol = None,
        loc_error_um2     = 0.0,
    ),
)


# =============================================================================
# AUTOLOAD CTCF BEDs AT IMPORT TIME
# =============================================================================
if AUTO_LOAD_CTCF_FROM_BED:
    for _cond_str, _rel in (("control",  CTCF_BED_CONTROL),
                             ("sorbitol", CTCF_BED_SORBITOL)):
        _abs = _os.path.join(_REPO_ROOT, _rel)
        if _os.path.exists(_abs):
            load_ctcf_from_bed(_abs, condition=_cond_str)
        else:
            import warnings as _warnings
            _warnings.warn(
                f"CTCF BED not found for condition='{_cond_str}': {_abs}\n"
                "Run setup_data.sh to generate it. Placeholder sites will be used.",
                UserWarning,
            )
