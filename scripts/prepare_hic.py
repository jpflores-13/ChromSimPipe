#!/usr/bin/env python
"""
Put one sample's Hi-C map into <outputDir>/mcool/<sample>.cool at the
config's `resolution`, whatever format it arrives in:

  .hic    hictk convert (.hic versions 8 and 9); every normalization in the
          file (KR, SCALE, VC, VC_SQRT, ...) becomes a bins column named
          after the method
  .mcool  the /resolutions/<resolution> cooler is copied out of the file
  .cool   copied as is, after checking its bin size is <resolution>

Everything downstream reads the single-resolution .cool this writes, so the
analysis does not depend on the input format. The normalization columns
found are printed, to check against analysis.hic_balance in the config.

Called by the convertHic rule in workflows/ChromSimPipe.snakefile:
    python scripts/prepare_hic.py --input <file> --output <out.cool> --resolution 1000
"""

import argparse
import json
import os
import shutil
import subprocess
import sys

FORMATS = (".hic", ".mcool", ".cool")


def hictk_bin():
    """The hictk executable: the one installed next to this Python (the
    simulation environment, even when it is not activated, as in
    check_inputs.sh), else the first on PATH."""
    local = os.path.join(os.path.dirname(sys.executable), "hictk")
    if os.access(local, os.X_OK):
        return local
    found = shutil.which("hictk")
    if found is None:
        raise FileNotFoundError("hictk is not installed in this environment "
                                "(it is in envs/cohesin_sim.yml; rebuild the environment)")
    return found


def run_hictk(*args, capture=True):
    """Run `hictk <args>` and return its stdout. HICTK_NO_TELEMETRY is set:
    hictk otherwise sends usage data to its developers' server on every call,
    which compute nodes without internet access would wait on."""
    env = {**os.environ, "HICTK_NO_TELEMETRY": "1"}
    cmd = [hictk_bin(), *map(str, args)]
    done = subprocess.run(cmd, env=env, capture_output=capture, text=True)
    if done.returncode != 0:
        detail = (done.stderr or "").strip().splitlines()
        raise RuntimeError(f"`hictk {' '.join(cmd[1:])}` failed"
                           + (f": {' '.join(l.strip() for l in detail[-3:])}" if detail else " (see the log above)"))
    return done.stdout


def hic_header(path):
    """
    What a .hic file's header says, read with hictk (format versions 8 and 9):
    {"version": int, "assembly": str, "resolutions": [bp, ...],
     "chromsizes": {name: length}}. Only the header is read, not the contacts.
    """
    meta = json.loads(run_hictk("metadata", "--output-format", "json", path))
    chroms = (line.split("\t") for line in run_hictk("dump", "--table", "chroms", path).splitlines()
              if line)
    return {
        "version": int(meta["format-version"]),
        "assembly": meta.get("assembly"),
        "resolutions": sorted(int(r) for r in meta["resolutions"]),
        "chromsizes": {name: int(length) for name, length in chroms},
    }


def hic_normalizations(path, resolution):
    """Normalization methods a .hic file stores at `resolution` (read from its footer)."""
    out = run_hictk("dump", "--table", "normalizations", "--resolution", resolution, path)
    return [line.strip() for line in out.splitlines() if line.strip()]


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
        header = hic_header(path)
        if resolution not in header["resolutions"]:
            sys.exit(f"{path} has no {resolution} bp resolution (it has: {header['resolutions']}). "
                     f"The simulation compares at `resolution` in the config ({resolution} bp, "
                     "one monomer per bin); use a .hic that includes it.")
        print(f"{path}: .hic version {header['version']}, assembly {header['assembly']}")
        ## Temporary files go next to the output, not to $TMPDIR (often small
        ## on compute nodes); --force replaces a partial file from a failed try.
        run_hictk("convert", path, out, "--resolutions", resolution,
                  "--tmpdir", os.path.dirname(os.path.abspath(out)), "--force", capture=False)
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
