"""
Simulation parameters loader - reads config/ChromSimConfig.yaml and the
samplesheet it names. Nothing here is hard-coded: every value comes from the
config (set CHROMSIM_CONFIG=/path/to/config.yaml to use another file).

Relative paths in the config and samplesheet are resolved against the repo
root, so one samplesheet pointing at data/... works on every cluster (each
cluster links its own files into data/). A relative outputDir is resolved
against $CHROMSIM_OUTPUT_ROOT when that is set (e.g. a scratch directory),
else against the repo root.

Several loci can be listed under `loci:`. Every value below describes ONE
active locus: $CHROMSIM_LOCUS if set (the workflow sets it for every job),
else the first locus in the config. LOCI and locus_paths() describe all of them.

Public API used by the Snakefile and scripts:

    CONFIG_PATH, REPO_ROOT, OUTPUT_ROOT, OUTPUT_DIR, PATHS, GENOME
    LOCI, LOCUS_NAMES, ACTIVE_LOCUS, locus_paths(), locus_tiling()
    CHROM, REGION_START, REGION_END, RESOLUTION, N_MONOMERS
    TILING, SIM_RUN, POLYSIM, N_REPLICATES, N_SHARDS
    SAMPLES, SIMULATION_CONDITIONS, ALL_PARAM_SETS
    GENOME_ASSEMBLY, CTCF_TYPES, CTCF_PEAKS, CTCF_BEDS, CTCF_SITES_BY_TYPE
    HIC_TYPES, HIC_FILES, MCOOL_FILES
    get_condition(), list_conditions()
    get_ctcf_arrays(), load_ctcf_from_bed()
    get_tiled_ctcf_arrays(), get_msd_pairs()
    ANALYSIS, DATA, MSD_PROBE, MSD_FIT, CALIBRATION

Sample names are built from the samplesheet columns listed in
``fileNamesFrom``; conditions refer to samples by that name.
"""

import csv as _csv
import numpy as np
import os as _os

# Locate and load config
REPO_ROOT   = _os.path.abspath(_os.path.join(_os.path.dirname(__file__), ".."))
CONFIG_PATH = _os.path.abspath(_os.path.expanduser(_os.environ.get(
    "CHROMSIM_CONFIG", _os.path.join(REPO_ROOT, "config", "ChromSimConfig.yaml"))))
_REPO_ROOT  = REPO_ROOT   # backwards-compatible alias

try:
    import yaml as _yaml
except ImportError:
    raise ImportError("PyYAML is required: pip install pyyaml (or conda install pyyaml)")

with open(CONFIG_PATH) as _f:
    _cfg = _yaml.safe_load(_f)


def _repo_path(path):
    """Expand ~ and resolve a config path relative to the repo root."""
    path = _os.path.expanduser(str(path))
    return path if _os.path.isabs(path) else _os.path.join(REPO_ROOT, path)


def _output_path(path):
    """Resolve outputDir: absolute as given, else under OUTPUT_ROOT."""
    path = _os.path.expanduser(str(path))
    return path if _os.path.isabs(path) else _os.path.join(OUTPUT_ROOT, path)


# Where a relative outputDir goes: $CHROMSIM_OUTPUT_ROOT (set per cluster,
# e.g. in ~/.bashrc, for a scratch filesystem) or the repo root.
OUTPUT_ROOT = _os.path.abspath(_os.path.expanduser(
    _os.environ.get("CHROMSIM_OUTPUT_ROOT") or REPO_ROOT))

# Genome FASTA (relative paths resolve against the repo root, e.g. data/...)
GENOME = _repo_path(_cfg["genome"])


# =============================================================================
# LOCI
# `loci:` is a list; each entry is simulated and analysed separately, in
# <outputDir>/loci/<name>/. Optional per locus: label and anchor_bp (MSD probe
# centre). An older single `locus:` block is still accepted.
# =============================================================================
RESOLUTION = int(_cfg.get("resolution", 1000))
_ASSEMBLY  = _cfg.get("genome_assembly", "hg38")
GENOME_ASSEMBLY = _ASSEMBLY

_loci_cfg = _cfg.get("loci")
if _loci_cfg is None and "locus" in _cfg:
    _loci_cfg = [_cfg["locus"]]
if not _loci_cfg:
    raise ValueError(f"No loci in {CONFIG_PATH}: add a `loci:` list.")

def _chrom_sizes(fasta):
    """Chromosome lengths from the genome's samtools index, if it exists yet
    (orientCTCF creates it on first use, before any GPU job)."""
    fai = f"{fasta}.fai"
    if not _os.path.isfile(fai):
        return {}
    with open(fai) as f:
        return {r[0]: int(r[1]) for r in (line.split("\t") for line in f) if len(r) > 1}


_CHROM_SIZES = _chrom_sizes(GENOME)

LOCI = {}
for _l in _loci_cfg:
    _name = str(_l["name"])
    if _name in LOCI or "/" in _name:
        raise ValueError(f"Locus names must be unique and contain no '/': {_name!r}")
    try:
        _start, _end = int(_l["region_start"]), int(_l["region_end"])
    except (TypeError, ValueError):
        raise ValueError(f"Locus {_name}: region_start and region_end must be whole numbers "
                         f"in bp (e.g. 64210000, not 64.21Mb or 64,210,000)") from None
    if _start < 0 or _end <= _start:
        raise ValueError(f"Locus {_name}: need 0 <= region_start < region_end (got {_start}-{_end})")
    if _start % RESOLUTION or _end % RESOLUTION:
        raise ValueError(
            f"Locus {_name}: region_start and region_end must be multiples of resolution "
            f"({RESOLUTION} bp) so each simulated monomer lines up with one Hi-C bin; "
            f"e.g. {_start // RESOLUTION * RESOLUTION}-{-(-_end // RESOLUTION) * RESOLUTION}")
    _chrom = str(_l["chrom"])
    if _CHROM_SIZES:
        if _chrom not in _CHROM_SIZES:
            raise ValueError(f"Locus {_name}: chromosome {_chrom!r} is not in {GENOME} "
                             f"(it has e.g. {list(_CHROM_SIZES)[:3]}); names must match exactly")
        if _end > _CHROM_SIZES[_chrom]:
            raise ValueError(f"Locus {_name}: region_end {_end} is past the end of {_chrom} "
                             f"({_CHROM_SIZES[_chrom]} bp); use at most "
                             f"{_CHROM_SIZES[_chrom] // RESOLUTION * RESOLUTION}")
    LOCI[_name] = dict(
        chrom        = str(_l["chrom"]),
        region_start = _start,
        region_end   = _end,
        region       = f"{_l['chrom']}:{_start}-{_end}",
        label        = _l.get("label", ""),
        anchor_bp    = int(_l.get("anchor_bp", (_start + _end) // 2)),
    )
LOCUS_NAMES = list(LOCI)

ACTIVE_LOCUS = _os.environ.get("CHROMSIM_LOCUS") or LOCUS_NAMES[0]
if ACTIVE_LOCUS not in LOCI:
    raise ValueError(f"CHROMSIM_LOCUS={ACTIVE_LOCUS!r} is not one of the config's loci {LOCUS_NAMES}")

_locus       = LOCI[ACTIVE_LOCUS]
CHROM        = _locus["chrom"]
REGION_START = _locus["region_start"]
REGION_END   = _locus["region_end"]
N_MONOMERS   = (REGION_END - REGION_START) // RESOLUTION
REGION       = _locus["region"]
ANCHOR_MONOMER = (_locus["anchor_bp"] - REGION_START) // RESOLUTION


# =============================================================================
# TILING
# Each locus is simulated as back-to-back copies on one long chain (more
# contact statistics per run; Yang et al. 2023 Nat Commun). Per locus:
#   tile_size = locus monomers + 2 * padding
#   n_tiles   = chain_size // tile_size, chrom_size = n_tiles * tile_size
# =============================================================================
_til = _cfg.get("tiling", {}) or {}
_PADDING    = int(_til.get("padding", 270))
_CHAIN_SIZE = int(_til.get("chain_size", 70000))


def locus_tiling(name):
    """Tiling of one locus (tile_size, n_tiles, chrom_size, padding)."""
    loc = LOCI[name]
    n = (loc["region_end"] - loc["region_start"]) // RESOLUTION
    tile = n + 2 * _PADDING
    n_tiles = _CHAIN_SIZE // tile
    if n_tiles < 1:
        raise ValueError(f"Locus {name} ({n} monomers + padding) is longer than "
                         f"tiling.chain_size ({_CHAIN_SIZE})")
    return dict(tile_size=tile, n_tiles=n_tiles, chrom_size=n_tiles * tile, padding=_PADDING)


TILING = locus_tiling(ACTIVE_LOCUS)


# =============================================================================
# SIMULATION RUN PARAMETERS
# =============================================================================
_run = _cfg.get("sim_run", {})
SIM_RUN = dict(
    seed_base            = int(_run.get("seed_base",            42)),
    hdf5_frames_per_file = int(_run.get("hdf5_frames_per_file", 100)),
    frames_per_replicate = int(_run.get("frames_per_replicate", 10000)),
    contact_radius       = float(_run.get("contact_radius",     3.0)),
    n_replicates         = int(_run.get("n_replicates", _cfg.get("n_replicates", 3))),
)
N_REPLICATES = SIM_RUN["n_replicates"]


# =============================================================================
# POLYSIM (the simulation model)
# Hansen lab polysim package (github.com/ahansenlab/polysim_bootcamp). Defaults
# are the bootcamp tutorial values (tutorial/sim3D_complete.py). A condition's
# `lifetime` is in LEF timesteps (= blocks), and each LEF leg steps with
# probability vlef per timestep, so the mean extruded length per leg is
# lifetime * vlef monomers.
# =============================================================================
_ps = _cfg.get("polysim", {}) or {}
POLYSIM = dict(
    density              = float(_ps.get("density",              0.3)),
    confinement_k        = float(_ps.get("confinement_k",        1.0)),
    repel                = float(_ps.get("repel",                3.0)),
    bond_wiggle          = float(_ps.get("bond_wiggle",          0.1)),
    angle_k              = float(_ps.get("angle_k",              0.05)),
    smc_bond_wiggle      = float(_ps.get("smc_bond_wiggle",      0.2)),
    smc_bond_dist        = float(_ps.get("smc_bond_dist",        0.5)),
    vlef                 = float(_ps.get("vlef",                 0.0025)),
    life_boost_stalled   = float(_ps.get("life_boost_stalled",   4.0)),
    platform             = str(_ps.get("platform",               "CUDA")),
    integrator           = str(_ps.get("integrator",             "langevin")),
    dt                   = int(_ps.get("dt",                     40)),
    colrate              = float(_ps.get("colrate",              0.01)),
    poly_steps_per_block = int(_ps.get("poly_steps_per_block",   33)),
    saveevery            = int(_ps.get("saveevery",              50)),
    initsteps            = int(_ps.get("initsteps",              540000)),
    initskip             = int(_ps.get("initskip",               0)),
    blocks_per_updater   = int(_ps.get("blocks_per_updater",     1000)),
)
N_SHARDS     = int(_cfg.get("n_shards", 2))


# =============================================================================
# OUTPUT LAYOUT (everything under outputDir)
# =============================================================================
OUTPUT_DIR = _output_path(_cfg.get("outputDir", "output"))


def locus_paths(name):
    """Output folders of one locus. Hi-C (.cool) and motifs are shared by all loci."""
    base = _os.path.join(OUTPUT_DIR, "loci", name)
    return dict(
        base       = base,
        mcool      = _os.path.join(OUTPUT_DIR, "mcool"),
        motifs     = _os.path.join(OUTPUT_DIR, "motifs"),
        ctcf       = _os.path.join(base, "ctcf"),
        sims       = _os.path.join(base, "sims"),
        analysis   = _os.path.join(base, "analysis"),
        figures    = _os.path.join(base, "figures"),
        done       = _os.path.join(base, "done"),
        logs       = _os.path.join(OUTPUT_DIR, "logs", name),
        benchmarks = _os.path.join(OUTPUT_DIR, "benchmarks"),
    )


PATHS = locus_paths(ACTIVE_LOCUS)


# =============================================================================
# SAMPLESHEET
# =============================================================================
SAMPLESHEET_PATH = _repo_path(_cfg["samplesheet"])
FILE_NAMES_FROM  = list(_cfg.get("fileNamesFrom") or [])
_MISSING = {"", "na", "none", "."}

with open(SAMPLESHEET_PATH, newline="") as _f:
    _rows = [
        {k.strip(): (v or "").strip() for k, v in r.items()}
        for r in _csv.DictReader(_f, delimiter="\t")
        if any((v or "").strip() for v in r.values())
    ]

for _col in ["CTCF_Peaks", "Hi-C"] + FILE_NAMES_FROM:
    if _rows and _col not in _rows[0]:
        raise KeyError(f"Samplesheet {SAMPLESHEET_PATH} has no '{_col}' column "
                       f"(columns: {list(_rows[0])})")
if not FILE_NAMES_FROM:
    raise ValueError("fileNamesFrom must list at least one samplesheet column")


def _value_or_none(v):
    return None if v.lower() in _MISSING else _repo_path(v)


SAMPLES = {}
for _r in _rows:
    _name = "_".join(_r[c] for c in FILE_NAMES_FROM)
    if _name in SAMPLES:
        raise ValueError(f"Duplicate sample name '{_name}' from fileNamesFrom "
                         f"{FILE_NAMES_FROM}; add a column that tells them apart.")
    SAMPLES[_name] = dict(
        metadata   = _r,
        ctcf_peaks = _value_or_none(_r["CTCF_Peaks"]),
        hic        = _value_or_none(_r["Hi-C"]),
    )

## Run name: metadata values shared by every sample, joined with '_'
## (the bagPipes `namer` pattern). Used to name run-level summary files.
_meta_cols = [c for c in (_rows[0] if _rows else {}) if c not in ("CTCF_Peaks", "Hi-C")]
_shared    = [next(iter(vals)) for vals in ({r[c] for r in _rows} for c in _meta_cols)
              if len(vals) == 1]
RUN_NAME   = "_".join(v for v in _shared if v) or "ChromSimPipe"


# =============================================================================
# CTCF BED FILE PATHS
# =============================================================================
# Lives in <outputDir>/loci/<locus>/ctcf/ctcf_oriented_{assembly}_{ctcf_type}.bed
LOCUS_TAG = f"{CHROM}_{REGION_START}_{REGION_END}"


def ctcf_bed_path(locus, ctcf_type):
    return _os.path.join(locus_paths(locus)["ctcf"],
                         f"ctcf_oriented_{_ASSEMBLY}_{ctcf_type}.bed")


def _ctcf_bed_path(ctcf_type):
    return ctcf_bed_path(ACTIVE_LOCUS, ctcf_type)

AUTO_LOAD_CTCF_FROM_BED = True

# Every ctcf_type declared by a condition, in first-seen order.
CTCF_TYPES = list(dict.fromkeys(
    c["ctcf_type"] for c in _cfg.get("conditions", []) if c.get("ctcf_type")
))
for _t in CTCF_TYPES:
    if _t not in SAMPLES or SAMPLES[_t]["ctcf_peaks"] is None:
        raise ValueError(
            f"Condition ctcf_type '{_t}' has no CTCF_Peaks in {SAMPLESHEET_PATH}. "
            f"Samples with peaks: {[k for k, v in SAMPLES.items() if v['ctcf_peaks']]}"
        )
CTCF_PEAKS = {t: SAMPLES[t]["ctcf_peaks"] for t in CTCF_TYPES}
CTCF_BEDS  = {t: _ctcf_bed_path(t) for t in CTCF_TYPES}


# =============================================================================
# CTCF SITE STORAGE (filled at import time by autoload below)
# =============================================================================
CTCF_SITES_BY_TYPE = {}


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
    )
    # optional per-condition overrides of the polysim: block
    for _k in ("vlef", "life_boost_stalled"):
        if _k in _c:
            _params[_k] = float(_c[_k])
    _cond = dict(
        name        = _c["name"],
        description = _c.get("description", _c["name"]),
        params      = _params,
        ctcf_type   = _c["ctcf_type"],
        compare_hic = _c.get("compare_hic", None),
    )
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
        hic = c["compare_hic"] or "-"
        print(
            f"{c['name']:<45} "
            f"{c['ctcf_type']:>10} "
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
        raise FileNotFoundError(
            f"No oriented CTCF BED loaded for ctcf_type='{condition}' "
            f"(expected {CTCF_BEDS.get(condition)}). Run the workflow's "
            "orientCTCF rule first."
        )
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


def get_tiled_ctcf_arrays(condition):
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
# Hi-C maps come from the samplesheet's Hi-C column; each is converted to
# <outputDir>/mcool/<sample>.cool.
HIC_FILES   = {k: v["hic"] for k, v in SAMPLES.items() if v["hic"]}
HIC_FORMATS = (".hic", ".mcool", ".cool")
for _k, _f in HIC_FILES.items():
    if not _f.lower().endswith(HIC_FORMATS):
        raise ValueError(f"Hi-C for sample '{_k}' is {_f}; the samplesheet's Hi-C column "
                         f"takes {', '.join(HIC_FORMATS)} files (or NA).")
HIC_TYPES   = list(HIC_FILES)
MCOOL_FILES = {t: _os.path.join(PATHS["mcool"], f"{t}.cool") for t in HIC_TYPES}

for _cond in SIMULATION_CONDITIONS:
    if _cond["compare_hic"] and _cond["compare_hic"] not in HIC_FILES:
        raise ValueError(
            f"Condition '{_cond['name']}' has compare_hic='{_cond['compare_hic']}', "
            f"but the samplesheet only has Hi-C maps for {HIC_TYPES}."
        )

DATA = dict(
    hic        = HIC_FILES,
    mcool      = MCOOL_FILES,
    resolution = RESOLUTION,
)


# =============================================================================
# MSD + CALIBRATION (analysis-side)
# =============================================================================
_ana = _cfg.get("analysis", {}) or {}
_apa = _ana.get("apa", {}) or {}
ANALYSIS = dict(
    n_jobs = int(_ana.get("n_jobs", 4)),
    msd    = bool(_ana.get("msd", False)),
    hic_balance = _ana.get("hic_balance", "auto"),
    apa    = dict(
        min_size_bins = int(_apa.get("min_size_bins", 20)),
        max_size_bins = int(_apa.get("max_size_bins", 1000)),
        pad           = int(_apa.get("pad", 10)),
    ),
)

_probe = _ana.get("msd_probe", {}) or {}
_auto  = _probe.get("auto", {}) or {}
MSD_PROBE = dict(
    # Explicit (monomer_a, monomer_b) pairs per ctcf_type; empty = auto-select.
    pairs  = dict(_probe.get("pairs") or {}),
    labels = dict(_probe.get("labels") or {}),
    auto = dict(
        anchor_bp              = _locus["anchor_bp"],
        min_sep_monomers       = int(_auto.get("min_sep_monomers", 200)),
        max_sep_monomers       = int(_auto.get("max_sep_monomers", 1500)),
        require_convergent     = bool(_auto.get("require_convergent", True)),
        n_pairs                = int(_auto.get("n_pairs", 2)),
        require_conserved      = bool(_auto.get("require_conserved", True)),
        conserved_tol_monomers = int(_auto.get("conserved_tol_monomers", 0)),
        prefer_flanking        = bool(_auto.get("prefer_flanking", True)),
        disjoint_anchors       = bool(_auto.get("disjoint_anchors", True)),
        relax_if_missing       = bool(_auto.get("relax_if_missing", True)),
    ),
)

_fit = _ana.get("msd_fit", {}) or {}
MSD_FIT = dict(
    lag_min_frames     = int(_fit.get("lag_min_frames", 1)),
    lag_max_frac       = float(_fit.get("lag_max_frac", 1 / 3)),
    fit_lag_min        = int(_fit.get("fit_lag_min", 5)),
    fit_lag_max_frac   = float(_fit.get("fit_lag_max_frac", 0.25)),
    min_n_lags_for_fit = int(_fit.get("min_n_lags_for_fit", 6)),
)

_cal     = _ana.get("calibration", {}) or {}
_cal_hic = _cal.get("hic", {}) or {}
_cal_msd = _cal.get("msd", {}) or {}
CALIBRATION = dict(
    default_anchor = _cal.get("default_anchor", "hic"),
    hic = dict(
        s_ref_bp             = int(_cal_hic.get("s_ref_bp", 100_000)),
        persistence_nm       = float(_cal_hic.get("persistence_nm", 50.0)),
        nm_per_bp            = float(_cal_hic.get("nm_per_bp", 0.34)),
        fixed_bp_per_monomer = RESOLUTION,
    ),
    msd = dict(
        expt_msd      = {k: _repo_path(v) for k, v in (_cal_msd.get("expt_msd") or {}).items()},
        loc_error_um2 = float(_cal_msd.get("loc_error_um2", 0.0)),
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
# AUTOLOAD CTCF BEDs AT IMPORT TIME
# =============================================================================
if AUTO_LOAD_CTCF_FROM_BED:
    # BEDs are workflow outputs, so they are absent when the Snakefile first
    # imports this module; get_ctcf_arrays() raises if one is still missing
    # when a script needs it.
    for _type, _path in CTCF_BEDS.items():
        if _os.path.exists(_path):
            load_ctcf_from_bed(_path, condition=_type)
