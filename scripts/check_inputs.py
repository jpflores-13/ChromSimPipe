#!/usr/bin/env python
"""
Check a run's inputs against each other before running anything.

For every sample in the samplesheet and every locus in the config:
  - the genome has each locus' chromosome, and the window fits inside it;
  - the CTCF peaks use the same chromosome names, have peaks in each window,
    and none lie past a chromosome's end (which would mean another assembly);
  - peak files named by a condition's sticky / loading / pausing options exist
    and have peaks in each window;
  - the Hi-C (.hic, .mcool or .cool) has the config's resolution, the same
    chromosome names, the same chromosome lengths as the genome (the same
    assembly), and, for .mcool/.cool, the normalization analysis.hic_balance
    asks for.

Prints a report and exits 1 if anything fails. Nothing is written.

    bash check_inputs.sh                      # finds the environment, runs this
    python scripts/check_inputs.py            # with the sim environment's Python
"""

import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

counts = defaultdict(int)


def report(status, msg):
    counts[status] += 1
    print(f"    [{status}] {msg}")


def other_style(name):
    """The same chromosome in the other naming style: chr1 <-> 1."""
    return name[3:] if name.startswith("chr") else f"chr{name}"


def name_hint(chrom, names):
    alt = other_style(chrom)
    return f" (it has {alt!r}: the files use different naming styles)" if alt in names else ""


def genome_sizes(fasta):
    """Chromosome lengths: from the .fai index if present, else by reading the FASTA."""
    fai = f"{fasta}.fai"
    if os.path.isfile(fai):
        with open(fai) as f:
            return {r[0]: int(r[1]) for r in (line.split("\t") for line in f) if len(r) > 1}
    if fasta.endswith(".gz"):
        raise ValueError("compressed FASTA without a .fai index; index it with samtools faidx")
    sizes, name = {}, None
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sizes[name] = 0
            elif name is not None:
                sizes[name] += len(line.strip())
    return sizes


def read_peaks(path):
    """Per chromosome: list of (start, end)."""
    peaks = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            r = line.split("\t")
            try:
                peaks[r[0]].append((int(r[1]), int(r[2])))
            except (IndexError, ValueError):
                continue
    return peaks


def hic_info(path, resolution):
    """(chromosome sizes, resolutions, assembly label, cooler for balance checks or None)."""
    name = path.lower()
    if name.endswith(".hic"):
        from hic2cool.hic2cool_utils import read_header
        with open(path, "rb") as f:
            chrs, resolutions, _, genome, _ = read_header(f)
        sizes = {c[1]: int(c[2]) for c in chrs.values() if str(c[1]).lower() != "all"}
        return sizes, sorted(resolutions), genome, None
    import cooler
    if name.endswith(".mcool"):
        groups = cooler.fileops.list_coolers(path)
        res = sorted(int(g.rsplit("/", 1)[1]) for g in groups
                     if g.startswith("/resolutions/") and g.rsplit("/", 1)[1].isdigit())
        uri = f"{path}::/resolutions/{resolution if resolution in res else res[0]}"
    else:
        uri, res = path, None
    clr = cooler.Cooler(uri)
    sizes = {k: int(v) for k, v in clr.chromsizes.items()}
    return sizes, res if res is not None else [clr.binsize], None, clr


def main():
    try:
        from configs import parameters as p
    except Exception as e:
        print(f"[FAIL] The config or samplesheet could not be read:\n  {e}")
        sys.exit(1)

    res = p.RESOLUTION
    print(f"Config: {p.CONFIG_PATH}")
    print(f"Resolution: {res} bp")
    print("Loci: " + ", ".join(f"{n} ({l['region']})" for n, l in p.LOCI.items()))

    print(f"\nGenome: {p.GENOME}")
    try:
        gsizes = genome_sizes(p.GENOME)
    except (OSError, ValueError) as e:
        report("FAIL", f"cannot read the genome: {e}")
        gsizes = {}
    if gsizes:
        print(f"    {len(gsizes)} sequences, named like {', '.join(list(gsizes)[:3])}")
    for name, loc in p.LOCI.items():
        c = loc["chrom"]
        if not gsizes:
            continue
        if c not in gsizes:
            report("FAIL", f"{name}: {c!r} is not in the genome{name_hint(c, gsizes)}")
        elif loc["region_end"] > gsizes[c]:
            report("FAIL", f"{name}: window ends at {loc['region_end']}, past the end of {c} ({gsizes[c]} bp)")
        else:
            report("PASS", f"{name}: {c} is {gsizes[c]:,} bp; the window fits")

    loci_chroms = {loc["chrom"] for loc in p.LOCI.values()}
    for sample, info in p.SAMPLES.items():
        print(f"\nSample {sample}")

        peaks_path = info["ctcf_peaks"]
        if peaks_path:
            print(f"  CTCF peaks: {peaks_path}")
            if not os.path.isfile(peaks_path):
                report("FAIL", "file not found")
            else:
                peaks = read_peaks(peaks_path)
                n = sum(len(v) for v in peaks.values())
                print(f"    {n:,} peaks on {len(peaks)} chromosomes, named like {', '.join(list(peaks)[:3])}")
                for name, loc in p.LOCI.items():
                    c = loc["chrom"]
                    if c not in peaks:
                        report("FAIL", f"{name}: no peaks on {c!r}{name_hint(c, peaks)}")
                        continue
                    inside = sum(1 for s, e in peaks[c] if e > loc["region_start"] and s < loc["region_end"])
                    report("PASS" if inside else "WARN",
                           f"{name}: {inside} peaks in the window" + ("" if inside else
                           " (no CTCF sites to simulate there)"))
                past = [c for c, v in peaks.items() if c in gsizes and max(e for _, e in v) > gsizes[c]]
                if past:
                    report("FAIL", f"peaks past the end of {', '.join(past[:3])} in the genome: "
                                   "the peaks and the genome come from different assemblies")
                elif gsizes:
                    report("PASS", "every peak lies within the genome's chromosomes")

        hic_path = info["hic"]
        if hic_path:
            print(f"  Hi-C: {hic_path}")
            if not os.path.isfile(hic_path):
                report("FAIL", "file not found")
                continue
            try:
                hsizes, hres, label, clr = hic_info(hic_path, res)
            except Exception as e:
                report("WARN", f"could not read it here ({e}); check again after convertHic, "
                               f"on {p.MCOOL_FILES[sample]}")
                continue
            if label:
                print(f"    assembly label in the file: {label}")
            report("PASS" if res in hres else "FAIL",
                   f"resolutions {hres}" + ("" if res in hres else f": {res} bp missing"))
            for c in sorted(loci_chroms):
                if c not in hsizes:
                    report("FAIL", f"{c!r} is not in the Hi-C{name_hint(c, hsizes)}")
                elif c in gsizes and hsizes[c] != gsizes[c]:
                    report("FAIL", f"{c} is {hsizes[c]:,} bp in the Hi-C but {gsizes[c]:,} bp in the "
                                   "genome: different assemblies")
                elif c in gsizes:
                    report("PASS", f"{c} has the same length in the Hi-C and the genome")
            if clr is not None and res in hres:
                from analysis.experimental_compare import resolve_balance
                try:
                    col = resolve_balance(clr, p.ANALYSIS["hic_balance"])
                    report("PASS", f"normalization: {col if col else 'raw counts'}")
                except ValueError as e:
                    report("FAIL", str(e))

    factor_files = {}
    for cond in p.SIMULATION_CONDITIONS:
        for label, path in p.condition_peak_files(cond).items():
            factor_files.setdefault(path, []).append(f"{cond['name']} {label}")
    for path, uses in factor_files.items():
        print(f"\nFactor peaks: {path}")
        print(f"    used by {', '.join(uses)}")
        if not os.path.isfile(path):
            report("FAIL", "file not found")
            continue
        peaks = read_peaks(path)
        print(f"    {sum(len(v) for v in peaks.values()):,} peaks, named like {', '.join(list(peaks)[:3])}")
        for name, loc in p.LOCI.items():
            c = loc["chrom"]
            if c not in peaks:
                report("FAIL", f"{name}: no peaks on {c!r}{name_hint(c, peaks)}")
                continue
            inside = sum(1 for st, e in peaks[c] if e > loc["region_start"] and st < loc["region_end"])
            report("PASS" if inside else "WARN", f"{name}: {inside} peaks in the window"
                   + ("" if inside else " (the option has no effect there)"))

    print(f"\n{counts['PASS']} passed, {counts['WARN']} warnings, {counts['FAIL']} failed")
    sys.exit(1 if counts["FAIL"] else 0)


if __name__ == "__main__":
    main()
