#!/usr/bin/env python
"""
Export an oriented BED6 file containing only the CTCF sites that are
conserved between mESC and cortical neurons (CN), using the
``_conserved_ctcf_sites`` helper already defined in ``configs/parameters.py``.

Why
---
The full per-cell-type CTCF BEDs (``CTCF_BED_MESC`` / ``CTCF_BED_NEURON`` in
``configs/parameters.py``) drive the simulation's CTCF barriers, the APA
pileup in ``analysis/absolute_quant.py``, and the contact-map CTCF arrow
overlays in ``analysis/ctcf_plotting.py``.  If instead of per-cell-type sites
you want every downstream analysis to operate on the mESC∩CN intersection
(same monomer index, same orientation), you need a single "conserved"
oriented BED that can be passed via ``--ctcf-bed-mesc`` *and*
``--ctcf-bed-neuron``.

This script builds exactly that BED in one pass.

How round-tripping works
------------------------
``_conserved_ctcf_sites`` returns ``(monomer_index, orientation)`` tuples.
We convert each monomer index back to genomic coordinates via

    summit_bp = REGION_START + idx * RESOLUTION + RESOLUTION // 2

and write a 1-bp BED record at ``[summit_bp, summit_bp + 1)``.  This
survives the ``(summit - REGION_START) // RESOLUTION`` round-trip used by
``parameters.load_ctcf_from_bed()``, so the exported file can be loaded
unchanged by the rest of the pipeline.

Usage
-----
    python scripts/export_conserved_ctcf_bed.py
    python scripts/export_conserved_ctcf_bed.py --tol 1       # allow 1-monomer drift
    python scripts/export_conserved_ctcf_bed.py \\
        --out data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed

Default output path:
    data/ctcf_oriented_CONSERVED_mESC_CN_chr3.bed
"""
from __future__ import annotations

import argparse
import os
import sys

# Make the repo root importable whether you call this from the repo root or
# from scripts/ directly.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from configs.parameters import (  # noqa: E402
    CHROM,
    REGION_START,
    RESOLUTION,
    N_MONOMERS,
    _conserved_ctcf_sites,
)


DEFAULT_OUT = os.path.join(
    "data", "ctcf_oriented_CONSERVED_mESC_CN_chr3.bed"
)


def build_bed_record(idx: int, orientation: int, name: str) -> str:
    """
    Build one BED6 line for a conserved CTCF site.

    The summit is placed at the centre of the monomer bin so that
    ``(summit - REGION_START) // RESOLUTION`` round-trips back to ``idx``.
    """
    summit_bp = REGION_START + idx * RESOLUTION + RESOLUTION // 2
    start_bp = summit_bp
    end_bp = summit_bp + 1
    strand = "+" if orientation > 0 else "-"
    # Score fixed at 600 (matches the example in README.md "Input file formats");
    # downstream code doesn't read the score column.
    return f"{CHROM}\t{start_bp}\t{end_bp}\t{name}\t600\t{strand}\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--tol", type=int, default=0,
        help="Monomer-index tolerance when matching mESC and neuron sites. "
             "0 = strict exact-match (default); 1-2 = accept sites that drift "
             "by one or two 1-kb bins between cell types.",
    )
    parser.add_argument(
        "--out", type=str, default=DEFAULT_OUT,
        help=f"Output BED path (default: {DEFAULT_OUT}).",
    )
    args = parser.parse_args()

    conserved = _conserved_ctcf_sites(tol_monomers=args.tol)
    if not conserved:
        print("ERROR: _conserved_ctcf_sites returned an empty list. "
              "Check that CTCF_BED_MESC and CTCF_BED_NEURON in "
              "configs/parameters.py are both loaded (AUTO_LOAD_CTCF_FROM_BED "
              "should be True).", file=sys.stderr)
        return 1

    # Sanity: all monomer indices should be inside the simulated region.
    bad = [i for i, _ in conserved if not (0 <= i < N_MONOMERS)]
    if bad:
        print(f"ERROR: {len(bad)} conserved sites have out-of-range monomer "
              f"indices (N_MONOMERS={N_MONOMERS}). Something is wrong with "
              f"the CTCF BEDs or REGION_START/RESOLUTION constants.",
              file=sys.stderr)
        return 1

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    plus = sum(1 for _, o in conserved if o > 0)
    minus = len(conserved) - plus

    with open(args.out, "w") as f:
        f.write("#chrom\tstart\tend\tname\tscore\tstrand\n")
        f.write(
            f"# Conserved oriented CTCF sites: mESC \u2229 CN "
            f"(tol={args.tol} monomer(s))\n"
        )
        f.write(
            f"# Source: configs.parameters._conserved_ctcf_sites() "
            f"intersecting CTCF_BED_MESC and CTCF_BED_NEURON\n"
        )
        f.write(
            f"# Region: {CHROM}:{REGION_START}-{REGION_START + N_MONOMERS * RESOLUTION}"
            f" at {RESOLUTION} bp/monomer\n"
        )
        f.write(
            f"# Counts: {len(conserved)} total "
            f"({plus} forward '+', {minus} reverse '-')\n"
        )
        for i, (idx, orient) in enumerate(conserved):
            name = f"CTCF_cons_{i:04d}"
            f.write(build_bed_record(idx, orient, name))

    print(
        f"Wrote {len(conserved)} conserved CTCF sites ({plus} +, {minus} -) "
        f"to {args.out}"
    )
    print(
        f"Region: {CHROM}:{REGION_START}-{REGION_START + N_MONOMERS * RESOLUTION}, "
        f"resolution: {RESOLUTION} bp/monomer, tolerance: {args.tol} monomer(s)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
