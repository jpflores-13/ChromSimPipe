"""
Simulation parameters loader — reads config/ChromSimConfig.yaml (falls back to SimConfig.yaml).

Presents the same public API as the original hardcoded parameters.py so all
downstream scripts work without modification:

    ACTIVE_LOCUS, LOCI, CHROM, REGION_START, REGION_END, RESOLUTION, N_MONOMERS
    POLYMER, SMC_BOND, TILING, SIM_RUN
    SIMULATION_CONDITIONS, ALL_PARAM_SETS
    GENOME_ASSEMBLY, CTCF_TYPES, CTCF_BEDS, CTCF_SITES_BY_TYPE
    HIC_TYPES, HIC_FILES, MCOOL_FILES
    PROMOTER_STALL_BED, PROMOTER_STALL_SITES
    get_condition(), list_conditions()
    get_ctcf_arrays(), load_ctcf_from_bed(), load_promoter_sites()
    get_tiled_ctcf_arrays(), get_msd_pairs()
    DATA, MSD_PROBE, MSD_FIT, CALIBRATION

CTCF types (e.g. "condition1", "condition2") and Hi-C types are whatever names the
config uses — nothing here assumes a particular pair of conditions.
"""

import numpy as np
import os as _os

# ── Locate and load config (ChromSimConfig.yaml preferred, falls back to SimConfig.yaml) ──
_REPO_ROOT   = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), ".."))
_CONFIG_PATH = _os.path.join(_REPO_ROOT, "config", "ChromSimConfig.yaml")
if not _os.path.exists(_CONFIG_PATH):
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
GENOME_ASSEMBLY = _ASSEMBLY
ANCHOR_MONOMER = (_locus["anchor_bp"] - REGION_START) // RESOLUTION


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

AUTO_LOAD_CTCF_FROM_BED = True

# Every ctcf_type declared by a condition, in first-seen order.
CTCF_TYPES = list(dict.fromkeys(
    c["ctcf_type"] for c in _cfg.get("conditions", []) if c.get("ctcf_type")
))
CTCF_BEDS = {t: _ctcf_bed_path(t) for t in CTCF_TYPES}

# Optional promoter stall-site BED, used by conditions with `ctcf_extra: promoter`.
# Relative paths resolve against the repo root.
PROMOTER_STALL_BED = _cfg.get("promoter_stall_bed") or None
if PROMOTER_STALL_BED and not _os.path.isabs(PROMOTER_STALL_BED):
    PROMOTER_STALL_BED = _os.path.join(_REPO_ROOT, PROMOTER_STALL_BED)


# =============================================================================
# CTCF SITE STORAGE (filled at import time by autoload below)
# =============================================================================
# Placeholder sites are used only for ctcf_types whose BED has not been
# generated yet (get_ctcf_arrays warns when that happens).
_PLACEHOLDER_SITES   = [(500, +1), (1000, -1), (1500, +1)]
CTCF_SITES_BY_TYPE   = {}
PROMOTER_STALL_SITES = []


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
        ctcf_type   = _c["ctcf_type"],
        compare_hic = _c.get("compare_hic", None),
    )
    if "ctcf_extra" in _c:
        if _c["ctcf_extra"] != "promoter":
            raise ValueError(
                f"Condition '{_c['name']}': unsupported ctcf_extra "
                f"'{_c['ctcf_extra']}' (only 'promoter' is implemented)."
            )
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

def get_ctcf_arrays(condition):
    """Return (positions, orientations) arrays for the given ctcf_type."""
    if condition not in CTCF_TYPES and condition not in CTCF_SITES_BY_TYPE:
        raise ValueError(
            f"Unknown ctcf_type '{condition}'. Types declared in the config: {CTCF_TYPES}"
        )
    sites = CTCF_SITES_BY_TYPE.get(condition)
    if sites is None:
        import warnings
        warnings.warn(
            f"Using PLACEHOLDER CTCF sites for ctcf_type='{condition}'! "
            "Run setup_data.sh to generate oriented BED files.",
            UserWarning,
        )
        sites = [s for s in _PLACEHOLDER_SITES if s[0] < N_MONOMERS]
    positions    = np.array([s[0] for s in sites], dtype=int)
    orientations = np.array([s[1] for s in sites], dtype=int)
    return positions, orientations


def load_ctcf_from_bed(bed_path, condition, region_start=None, resolution=None):
    """
    Load oriented CTCF sites from a BED6 file into CTCF_SITES_BY_TYPE[condition].
    """
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
    CTCF_SITES_BY_TYPE[condition] = sites
    return sites


def load_promoter_sites(bed_path, region_start=None, resolution=None):
    """
    Load promoter-based cohesin stall sites.

    Promoter anchors are non-directional: each site gets orientation 0, which
    LEFSimulator treats as a barrier to both extrusion directions.
    """
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
            if 0 <= monomer < N_MONOMERS:
                sites.append((monomer, 0))
    sites = sorted(set(sites))
    PROMOTER_STALL_SITES = sites
    return sites


def get_barrier_arrays(condition, ctcf_extra=None):
    """
    Return (positions, orientations) for one ctcf_type, plus any extra
    barriers the condition asks for (currently only ctcf_extra="promoter").
    """
    positions, orientations = get_ctcf_arrays(condition)
    if ctcf_extra == "promoter":
        if not PROMOTER_STALL_SITES:
            raise RuntimeError(
                "Condition requests ctcf_extra='promoter' but no promoter stall "
                f"sites are loaded. Set promoter_stall_bed in {_CONFIG_PATH} "
                f"(currently: {PROMOTER_STALL_BED!r}) to a BED with sites in "
                f"{CHROM}:{REGION_START}-{REGION_END}."
            )
        positions = np.concatenate(
            [positions, np.array([s[0] for s in PROMOTER_STALL_SITES], dtype=int)])
        orientations = np.concatenate(
            [orientations, np.array([s[1] for s in PROMOTER_STALL_SITES], dtype=int)])
    elif ctcf_extra is not None:
        raise ValueError(f"Unsupported ctcf_extra '{ctcf_extra}'")
    return positions, orientations


def get_tiled_ctcf_arrays(condition, ctcf_extra=None):
    """Return barrier positions/orientations tiled across the full simulated chromosome."""
    positions, orientations = get_barrier_arrays(condition, ctcf_extra)
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
# Hi-C maps come from the `hic:` mapping in the config ({hic_type: path}).
# Older configs with top-level `hic_<type>:` keys are still accepted.
HIC_FILES = dict(_cfg.get("hic") or {})
for _k, _v in _cfg.items():
    if _k.startswith("hic_") and _k[4:] not in HIC_FILES:
        HIC_FILES[_k[4:]] = _v
HIC_TYPES   = list(HIC_FILES)
MCOOL_FILES = {t: f"data/mcool/{t}.cool" for t in HIC_TYPES}

for _cond in SIMULATION_CONDITIONS:
    if _cond["compare_hic"] and _cond["compare_hic"] not in HIC_FILES:
        raise ValueError(
            f"Condition '{_cond['name']}' has compare_hic='{_cond['compare_hic']}', "
            f"but the config only defines Hi-C maps for {HIC_TYPES}."
        )

DATA = dict(
    hic        = HIC_FILES,
    mcool      = MCOOL_FILES,
    resolution = RESOLUTION,
)


# =============================================================================
# MSD + CALIBRATION (analysis-side)
# =============================================================================
MSD_PROBE = dict(
    # Explicit (monomer_a, monomer_b) pairs per ctcf_type; empty = auto-select.
    pairs  = {},
    labels = {},
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
        expt_msd      = {},   # {ctcf_type: CSV path}
        loc_error_um2 = 0.0,
    ),
)


def get_msd_pairs(ctcf_type):
    """
    Return (pairs, labels) of CTCF anchor monomers to probe for two-point MSD.

    Uses MSD_PROBE["pairs"][ctcf_type] if set; otherwise auto-selects
    convergent CTCF pairs around the locus anchor (see MSD_PROBE["auto"]).
    """
    explicit = MSD_PROBE["pairs"].get(ctcf_type)
    if explicit:
        pairs = [tuple(p) for p in explicit]
        labels = MSD_PROBE["labels"].get(ctcf_type)
        if not labels or len(labels) != len(pairs):
            labels = [f"{a}-{b}" for a, b in pairs]
        return pairs, list(labels)
    pairs = _auto_select_msd_pairs(ctcf_type)
    return pairs, [f"auto:{a}-{b}" for a, b in pairs]


def _conserved_ctcf_sites(tol_monomers=0):
    """
    Return (monomer, orientation) CTCF sites present in EVERY ctcf_type,
    matched within tol_monomers. The first type's monomer index is kept.
    """
    import bisect
    per_type = []
    for t in CTCF_TYPES:
        pos, ori = get_ctcf_arrays(t)
        per_type.append(list(zip(pos.tolist(), ori.tolist())))
    if not per_type:
        return []
    conserved = sorted(set(per_type[0]))
    for other in per_type[1:]:
        if tol_monomers <= 0:
            conserved = sorted(set(conserved) & set(other))
            continue
        by_ori = {}
        for p, o in other:
            by_ori.setdefault(o, []).append(p)
        for lst in by_ori.values():
            lst.sort()
        kept = []
        for p, o in conserved:
            lst = by_ori.get(o, [])
            i = bisect.bisect_left(lst, p)
            near = lst[max(i - 1, 0):i + 1]
            if any(abs(p - c) <= tol_monomers for c in near):
                kept.append((p, o))
        conserved = kept
    return conserved


def _auto_select_msd_pairs(ctcf_type):
    """
    Pick MSD_PROBE["auto"]["n_pairs"] CTCF pairs separated by
    [min_sep, max_sep] monomers, preferring convergent (+, -) pairs that
    flank the locus anchor and, if require_conserved, sites shared by all
    ctcf_types. Constraints relax in that order if nothing matches.
    """
    auto = MSD_PROBE["auto"]
    anchor_monomer = (auto["anchor_bp"] - REGION_START) // RESOLUTION
    min_sep, max_sep = auto["min_sep_monomers"], auto["max_sep_monomers"]
    n_pairs = auto["n_pairs"]
    prefer_flanking = auto.get("prefer_flanking", True)
    disjoint = auto.get("disjoint_anchors", True)
    require_conserved = auto.get("require_conserved", True)
    tol = auto.get("conserved_tol_monomers", 0)

    def _build_sites(require_cons):
        if require_cons:
            return _conserved_ctcf_sites(tol_monomers=tol)
        pos, ori = get_ctcf_arrays(ctcf_type)
        return list(zip(pos.tolist(), ori.tolist()))

    def _enumerate(sites, require_conv):
        cands = []
        for a_idx, (pa, oa) in enumerate(sites):
            for pb, ob in sites[a_idx + 1:]:
                sep = pb - pa
                if sep < min_sep or sep > max_sep:
                    continue
                if require_conv and (oa != +1 or ob != -1):
                    continue
                flanks = pa < anchor_monomer < pb
                rank = (0 if (flanks and prefer_flanking) else 1,
                        abs((pa + pb) / 2.0 - anchor_monomer))
                cands.append((rank, pa, pb))
        cands.sort(key=lambda r: r[0])
        return cands

    sites = _build_sites(require_conserved)
    candidates = _enumerate(sites, auto["require_convergent"])
    if not candidates and require_conserved and auto["relax_if_missing"]:
        sites = _build_sites(False)
        candidates = _enumerate(sites, auto["require_convergent"])
    if not candidates and auto["relax_if_missing"]:
        candidates = _enumerate(sites, False)
    if not candidates:
        raise RuntimeError(
            f"Could not auto-select an MSD pair for ctcf_type '{ctcf_type}': no CTCF "
            f"pair satisfies min_sep={min_sep}, max_sep={max_sep}. Set explicit "
            f"pairs in MSD_PROBE['pairs']['{ctcf_type}']."
        )

    chosen, used = [], set()
    for _, pa, pb in candidates:
        if disjoint and (pa in used or pb in used):
            continue
        chosen.append((int(pa), int(pb)))
        used.update((pa, pb))
        if len(chosen) >= n_pairs:
            break
    for _, pa, pb in candidates:
        if len(chosen) >= n_pairs:
            break
        if (int(pa), int(pb)) not in chosen:
            chosen.append((int(pa), int(pb)))
    return chosen


# =============================================================================
# AUTOLOAD CTCF + PROMOTER BEDs AT IMPORT TIME
# =============================================================================
if AUTO_LOAD_CTCF_FROM_BED:
    import warnings as _warnings
    for _type, _rel in CTCF_BEDS.items():
        _abs = _os.path.join(_REPO_ROOT, _rel)
        if _os.path.exists(_abs):
            load_ctcf_from_bed(_abs, condition=_type)
        else:
            _warnings.warn(
                f"CTCF BED not found for ctcf_type='{_type}': {_abs}\n"
                "Run setup_data.sh to generate it. Placeholder sites will be used.",
                UserWarning,
            )
    if PROMOTER_STALL_BED and _os.path.exists(PROMOTER_STALL_BED):
        load_promoter_sites(PROMOTER_STALL_BED)
    elif any(c.get("ctcf_extra") == "promoter" for c in SIMULATION_CONDITIONS):
        _warnings.warn(
            "Some conditions use ctcf_extra='promoter' but promoter_stall_bed "
            f"is not set or missing ({PROMOTER_STALL_BED!r}). Those conditions "
            "are skipped by Snakemake and refuse to simulate until it is provided.",
            UserWarning,
        )
