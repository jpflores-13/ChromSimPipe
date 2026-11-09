#!/usr/bin/env python
"""
Put one sample's Hi-C map into <outputDir>/mcool/<sample>.cool at the
config's `resolution`, whatever format it arrives in:

  .hic    hic2cool convert; Juicer's normalizations become the KR, VC and
          VC_SQRT columns of the bins table
  .mcool  the /resolutions/<resolution> cooler is copied out of the file
  .cool   copied as is, after checking its bin size is <resolution>

Everything downstream reads the single-resolution .cool this writes, so the
analysis does not depend on the input format. The normalization columns
found are printed, to check against analysis.hic_balance in the config.

Called by the convertHic rule in workflows/ChromSimPipe.snakefile:
    python scripts/prepare_hic.py --input <file> --output <out.cool> --resolution 1000
"""

import argparse
import os
import subprocess
import sys

FORMATS = (".hic", ".mcool", ".cool")


def mcool_resolutions(path):
    """Bin sizes (bp) stored under /resolutions/ in an .mcool file."""
    import cooler
    out = []
    for group in cooler.fileops.list_coolers(path):
        parts = group.strip("/").split("/")
        if len(parts) == 2 and parts[0] == "resolutions" and parts[1].isdigit():
            out.append(int(parts[1]))
    return sorted(out)


def prepare(path, out, resolution):
    name = path.lower()
    if not name.endswith(FORMATS):
        sys.exit(f"{path}: the samplesheet's Hi-C column takes {', '.join(FORMATS)} files")
    if not os.path.exists(path):
        sys.exit(f"{path} does not exist")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    if name.endswith(".hic"):
        subprocess.run(["hic2cool", "convert", path, out, "-r", str(resolution)], check=True)
    else:
        import cooler
        if name.endswith(".mcool"):
            have = mcool_resolutions(path)
            if resolution not in have:
                sys.exit(f"{path} has no {resolution} bp resolution (it has: {have or 'none'}). "
                         f"The simulation compares at `resolution` in the config ({resolution} bp, "
                         "one monomer per bin); use an .mcool that includes it.")
            src = f"{path}::/resolutions/{resolution}"
        else:
            src = path
            binsize = cooler.Cooler(src).binsize
            if binsize != resolution:
                sys.exit(f"{path} is binned at {binsize} bp, but `resolution` in the config is "
                         f"{resolution} bp; give a .cool (or .mcool) at {resolution} bp.")
        cooler.fileops.cp(src, out, overwrite=True)

    import cooler
    clr = cooler.Cooler(out)
    columns = [c for c in clr.bins().columns if c not in ("chrom", "start", "end")]
    print(f"{out}: {clr.binsize} bp bins, normalizations available: {columns or 'none'}")


def main():
    ap = argparse.ArgumentParser(description="Convert or copy one Hi-C map to a single-resolution .cool")
    ap.add_argument("--input", required=True, help=".hic, .mcool or .cool file")
    ap.add_argument("--output", required=True, help="Output .cool path")
    ap.add_argument("--resolution", type=int, required=True, help="Bin size in bp")
    args = ap.parse_args()
    prepare(args.input, args.output, args.resolution)


if __name__ == "__main__":
    main()
