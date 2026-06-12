#!/usr/bin/env python3
"""
Validate CTCF BED files: parsing, coordinate conversion, and polychrom format.

BED paths and conditions are read from configs/parameters.py (which reads
ChromSimConfig.yaml). No arguments needed for a standard run.

Usage:
    python scripts/validate_ctcf.py
"""

import os
import sys
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from configs.parameters import (
    CHROM, REGION_START, REGION_END, RESOLUTION, N_MONOMERS,
    ANCHOR_MONOMER, CTCF_BEDS,
    CTCF_SITES_PLACEHOLDER,
    get_ctcf_arrays, load_ctcf_from_bed,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
PASS = "\033[92m PASS \033[0m"
FAIL = "\033[91m FAIL \033[0m"
WARN = "\033[93m WARN \033[0m"
INFO = "\033[94m INFO \033[0m"

n_pass = 0
n_fail = 0
n_warn = 0


def check(condition, msg, warn_only=False):
    global n_pass, n_fail, n_warn
    if condition:
        print(f"  [{PASS}] {msg}")
        n_pass += 1
        return True
    elif warn_only:
        print(f"  [{WARN}] {msg}")
        n_warn += 1
        return False
    else:
        print(f"  [{FAIL}] {msg}")
        n_fail += 1
        return False


def info(msg):
    print(f"  [{INFO}] {msg}")


def header(title):
    print()
    print(f"{'─' * 70}")
    print(f"  {title}")
    print(f"{'─' * 70}")


# ---------------------------------------------------------------------------
# 1. Validate a single BED file (raw parsing)
# ---------------------------------------------------------------------------
def validate_bed_file(bed_path, label):
    header(f"BED FILE: {label} — {os.path.basename(bed_path)}")

    if not check(os.path.isfile(bed_path), f"File exists: {bed_path}"):
        return None

    with open(bed_path) as f:
        raw_lines = f.readlines()

    info(f"Total lines: {len(raw_lines)}")

    has_track = any(l.strip().startswith("track") for l in raw_lines)
    if has_track:
        info("Track header found (will be skipped by loader)")

    data_lines = []
    bad_lines = []
    for i, line in enumerate(raw_lines, 1):
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("track"):
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            bad_lines.append((i, line, "fewer than 6 columns"))
            continue
        try:
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[5]
            data_lines.append({
                "line": i, "chrom": chrom, "start": start, "end": end,
                "strand": strand, "raw": line,
            })
        except (ValueError, IndexError) as e:
            bad_lines.append((i, line, str(e)))

    check(len(bad_lines) == 0,
          f"All data lines parse correctly ({len(data_lines)} sites, {len(bad_lines)} bad)")
    if bad_lines:
        for ln, content, reason in bad_lines[:5]:
            print(f"        Line {ln}: {reason}")
            print(f"        → {content[:80]}")

    if not data_lines:
        check(False, "At least one data line found")
        return None

    on_chrom = [d for d in data_lines if d["chrom"] == CHROM]
    off_chrom = [d for d in data_lines if d["chrom"] != CHROM]
    check(len(on_chrom) > 0, f"Sites on {CHROM}: {len(on_chrom)}")
    if off_chrom:
        check(False, f"Sites on other chromosomes: {len(off_chrom)} (will be ignored)",
              warn_only=True)

    strands = set(d["strand"] for d in on_chrom)
    check(strands <= {"+", "-"}, f"Strand values are +/- only: {strands}")
    n_plus = sum(1 for d in on_chrom if d["strand"] == "+")
    n_minus = sum(1 for d in on_chrom if d["strand"] == "-")
    info(f"Orientation breakdown: {n_plus} forward (+), {n_minus} reverse (-)")

    in_region = []
    out_of_region = []
    for d in on_chrom:
        summit = (d["start"] + d["end"]) // 2
        if REGION_START <= summit < REGION_END:
            in_region.append(d)
        else:
            out_of_region.append(d)

    check(len(in_region) > 0,
          f"Sites within simulation region [{REGION_START:,}-{REGION_END:,}]: {len(in_region)}")
    if out_of_region:
        info(f"Sites outside region (will be ignored): {len(out_of_region)}")

    bad_coords = [d for d in on_chrom if d["start"] >= d["end"]]
    check(len(bad_coords) == 0,
          f"All peaks have start < end: {len(bad_coords)} violations")

    widths = [d["end"] - d["start"] for d in on_chrom]
    if widths:
        info(f"Peak widths: min={min(widths)} bp, max={max(widths)} bp, "
             f"median={sorted(widths)[len(widths)//2]} bp")
        check(all(10 < w < 10000 for w in widths),
              "Peak widths are reasonable (10-10,000 bp)", warn_only=True)

    return in_region


# ---------------------------------------------------------------------------
# 2. Validate conversion to monomer coordinates
# ---------------------------------------------------------------------------
def validate_conversion(sites_data, label):
    header(f"COORDINATE CONVERSION: {label}")

    if sites_data is None:
        check(False, "No sites to convert (BED parsing failed)")
        return None

    monomers = []
    for d in sites_data:
        summit = (d["start"] + d["end"]) // 2
        monomer = (summit - REGION_START) // RESOLUTION
        ori = +1 if d["strand"] == "+" else -1
        monomers.append((monomer, ori))

    monomers.sort(key=lambda x: x[0])

    positions = [m[0] for m in monomers]
    orientations = [m[1] for m in monomers]

    info(f"Converted {len(monomers)} sites to monomer coordinates")

    check(all(0 <= p < N_MONOMERS for p in positions),
          f"All positions in [0, {N_MONOMERS}): "
          f"range [{min(positions)}, {max(positions)}]")

    unique_pos = set(positions)
    if len(unique_pos) < len(positions):
        dups = len(positions) - len(unique_pos)
        check(False, f"No duplicate monomer positions: {dups} duplicates found",
              warn_only=True)
        info("(Multiple peaks mapping to the same monomer — only last one used)")

    check(all(o in (+1, -1) for o in orientations),
          f"All orientations are +1 or -1")
    check(positions == sorted(positions),
          "Positions are sorted in ascending order")

    info(f"Monomer positions: {positions}")
    info(f"Orientations:      {['+' if o == 1 else '-' for o in orientations]}")

    density = len(monomers) / (N_MONOMERS / 100)
    info(f"CTCF density: {density:.1f} sites per 100 monomers (100 kb)")

    convergent = 0
    for i in range(len(monomers) - 1):
        if monomers[i][1] == +1 and monomers[i + 1][1] == -1:
            convergent += 1
    info(f"Convergent pairs (+1, -1): {convergent} potential TAD boundaries")

    return monomers


# ---------------------------------------------------------------------------
# 3. Validate load_ctcf_from_bed() function
# ---------------------------------------------------------------------------
def validate_loader(bed_path, condition, expected_manual):
    header(f"LOADER FUNCTION: load_ctcf_from_bed(condition='{condition}')")

    try:
        loaded = load_ctcf_from_bed(bed_path, condition=condition)
        check(True, "load_ctcf_from_bed() ran without errors")
    except Exception as e:
        check(False, f"load_ctcf_from_bed() raised: {e}")
        return

    check(isinstance(loaded, list), f"Returns a list: {type(loaded).__name__}")
    check(len(loaded) > 0, f"Returned {len(loaded)} sites")

    if expected_manual is not None:
        check(len(loaded) == len(expected_manual),
              f"Site count matches manual parse: {len(loaded)} vs {len(expected_manual)}")

        if len(loaded) == len(expected_manual):
            match = all(a == b for a, b in zip(loaded, expected_manual))
            check(match, "Sites match manual conversion exactly")
            if not match:
                for i, (a, b) in enumerate(zip(loaded, expected_manual)):
                    if a != b:
                        print(f"        Mismatch at index {i}: loader={a} vs manual={b}")

    if loaded:
        check(all(isinstance(s, tuple) and len(s) == 2 for s in loaded),
              "All entries are (position, orientation) tuples")
        check(all(isinstance(s[0], int) and isinstance(s[1], int) for s in loaded),
              "Positions and orientations are integers")

    return loaded


# ---------------------------------------------------------------------------
# 4. Validate get_ctcf_arrays() output
# ---------------------------------------------------------------------------
def validate_arrays(condition, expected_sites):
    header(f"POLYCHROM FORMAT: get_ctcf_arrays(condition='{condition}')")

    try:
        positions, orientations = get_ctcf_arrays(condition=condition)
        check(True, "get_ctcf_arrays() ran without errors")
    except Exception as e:
        check(False, f"get_ctcf_arrays() raised: {e}")
        return

    check(isinstance(positions, np.ndarray), f"Positions: numpy array ({positions.dtype})")
    check(isinstance(orientations, np.ndarray), f"Orientations: numpy array ({orientations.dtype})")
    check(positions.dtype in (np.int32, np.int64, int), f"Position dtype is int: {positions.dtype}")
    check(orientations.dtype in (np.int32, np.int64, int), f"Orientation dtype is int: {orientations.dtype}")
    check(len(positions) == len(orientations),
          f"Arrays have same length: {len(positions)} positions, {len(orientations)} orientations")
    check(np.all(positions >= 0) and np.all(positions < N_MONOMERS),
          f"All positions in [0, {N_MONOMERS})")
    check(set(orientations.tolist()) <= {+1, -1},
          f"Orientations are +1/-1 only: {set(orientations.tolist())}")

    if expected_sites:
        expected_pos = np.array([s[0] for s in expected_sites], dtype=int)
        expected_ori = np.array([s[1] for s in expected_sites], dtype=int)
        check(np.array_equal(positions, expected_pos), "Positions match loaded sites")
        check(np.array_equal(orientations, expected_ori), "Orientations match loaded sites")

    left_stall = np.zeros(N_MONOMERS, dtype=bool)
    right_stall = np.zeros(N_MONOMERS, dtype=bool)
    for pos, ori in zip(positions, orientations):
        if ori == +1:
            left_stall[pos] = True
        elif ori == -1:
            right_stall[pos] = True

    info(f"LEF stall arrays: {left_stall.sum()} left-blockers (+1), {right_stall.sum()} right-blockers (-1)")
    info(f"Total occupied monomers: {(left_stall | right_stall).sum()} / {N_MONOMERS}")

    # Check CTCF near the locus anchor point
    anchor_range = range(ANCHOR_MONOMER - 50, ANCHOR_MONOMER + 50)
    near_anchor = [(p, o) for p, o in zip(positions, orientations) if p in anchor_range]
    if near_anchor:
        anchor_bp = ANCHOR_MONOMER * RESOLUTION + REGION_START
        info(f"CTCF sites within ±50 kb of anchor (monomer {ANCHOR_MONOMER}, {anchor_bp:,} bp):")
        for p, o in near_anchor:
            dist = p - ANCHOR_MONOMER
            strand = "+" if o == +1 else "-"
            genomic = p * RESOLUTION + REGION_START
            info(f"  monomer {p} ({genomic:,} bp), orientation {strand}, "
                 f"{'upstream' if dist < 0 else 'downstream'} {abs(dist)} kb")
    else:
        check(False, f"No CTCF sites near anchor ±50 kb (monomer {ANCHOR_MONOMER})", warn_only=True)


# ---------------------------------------------------------------------------
# 5. Compare two conditions
# ---------------------------------------------------------------------------
def compare_conditions(cond_a, cond_b):
    header(f"COMPARISON: {cond_a} vs {cond_b} CTCF landscapes")

    pos_a, ori_a = get_ctcf_arrays(condition=cond_a)
    pos_b, ori_b = get_ctcf_arrays(condition=cond_b)

    info(f"{cond_a}: {len(pos_a)} sites")
    info(f"{cond_b}: {len(pos_b)} sites")

    tolerance = 2  # monomers = 2 kb
    shared = 0
    a_only = 0
    b_only = 0

    for p in pos_a:
        if any(abs(p - q) <= tolerance for q in pos_b):
            shared += 1
        else:
            a_only += 1

    for p in pos_b:
        if not any(abs(p - q) <= tolerance for q in pos_a):
            b_only += 1

    info(f"Shared (within {tolerance} kb): {shared}")
    info(f"{cond_a}-specific: {a_only}")
    info(f"{cond_b}-specific: {b_only}")

    check(not np.array_equal(pos_a, pos_b),
          f"{cond_a} and {cond_b} have DIFFERENT CTCF landscapes")

    anchor_window = range(ANCHOR_MONOMER - 100, ANCHOR_MONOMER + 100)
    a_near = sum(1 for p in pos_a if p in anchor_window)
    b_near = sum(1 for p in pos_b if p in anchor_window)
    info(f"Sites within ±100 kb of anchor: {cond_a}={a_near}, {cond_b}={b_near}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("  CTCF VALIDATION — Pre-simulation checks")
    print("=" * 70)
    print(f"  Region:     {CHROM}:{REGION_START:,}-{REGION_END:,}")
    print(f"  Resolution: {RESOLUTION} bp/monomer")
    print(f"  Monomers:   {N_MONOMERS}")
    print(f"  Anchor:     monomer {ANCHOR_MONOMER} ({ANCHOR_MONOMER * RESOLUTION + REGION_START:,} bp)")
    print(f"  Conditions: {list(CTCF_BEDS)}")
    print(f"  Placeholder CTCF: {CTCF_SITES_PLACEHOLDER}")

    loaded_per_condition = {}

    for condition, bed_path in CTCF_BEDS.items():
        if os.path.isfile(bed_path):
            raw = validate_bed_file(bed_path, condition)
            monomers = validate_conversion(raw, condition)
            loaded = validate_loader(bed_path, condition, monomers)
            validate_arrays(condition, loaded)
            loaded_per_condition[condition] = loaded
        else:
            header(f"{condition.upper()} — BED FILE NOT FOUND")
            info(f"Expected: {bed_path}")
            info("Using placeholder sites from parameters.py")
            validate_arrays(condition, None)

    # Compare the first two conditions if both were loaded
    conditions = list(CTCF_BEDS)
    if len(conditions) >= 2:
        cond_a, cond_b = conditions[0], conditions[1]
        if cond_a in loaded_per_condition and cond_b in loaded_per_condition:
            compare_conditions(cond_a, cond_b)

    header("SUMMARY")
    print(f"  Passed: {n_pass}")
    print(f"  Warnings: {n_warn}")
    print(f"  Failed: {n_fail}")
    print()

    if n_fail > 0:
        print("  \033[91m✗ CTCF validation FAILED — fix errors before running simulations\033[0m")
        sys.exit(1)
    elif n_warn > 0:
        print("  \033[93m⚠ CTCF validation passed with warnings — review before production\033[0m")
    else:
        print("  \033[92m✓ All checks passed — CTCF sites are ready for simulation\033[0m")
    print()


if __name__ == "__main__":
    main()
