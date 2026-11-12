#!/usr/bin/env python
"""
CTCF occupancy per oriented site, from the ChIP bigWig and, optionally, an input
bigWig (ctcfOccupancy rule).

For every CTCF peak on the locus chromosome:
  S = mean ChIP signal over summit +/- window/2
  I = mean input signal over summit +/- input_window/2 (local background)
  E = max(S - I, 0)                      enrichment above input (I = 0 without an input)
  occupancy = clip(E / ref, 0, 1)        ref = the chromosome's max E (or a percentile)
so the strongest peak on the chromosome is 1 and a site at input level is 0
(Michele's scaling). Each oriented site in the locus BED gets the occupancy of the
peak it came from (matched by peak name). The simulation turns occupancy into a
stall probability (ctcf_capture {by: occupancy}); several sites on one monomer
combine as 1 - prod(1 - p) (Harvey's supplement, section 8.3).

Output TSV, one row per oriented site, in BED order:
  name chrom start end strand chip input enrichment occupancy
plus header comments with the chromosome-wide reference.
"""
import argparse, math
import numpy as np
import pyBigWig


def mean_signal(bw, chrom, a, b):
    a = max(0, a)
    b = min(b, bw.chroms(chrom))
    if b <= a:
        return 0.0
    v = bw.stats(chrom, a, b, type="mean")[0]
    return 0.0 if v is None or math.isnan(v) else float(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peaks", required=True, help="CTCF narrowPeak (whole genome or chromosome)")
    ap.add_argument("--sites", required=True, help="oriented CTCF BED for the locus (orientCTCF)")
    ap.add_argument("--chip", required=True, help="CTCF ChIP bigWig")
    ap.add_argument("--input", nargs="?", default=None, const=None,
                    help="input bigWig (optional; passed empty or omitted = ChIP signal alone)")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--window", type=int, default=200, help="ChIP window around the summit (bp)")
    ap.add_argument("--input-window", type=int, default=2000, help="input window around the summit (bp)")
    ap.add_argument("--ref", default="max", help="'max' or a percentile like 'p99' of chromosome enrichment")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    chip = pyBigWig.open(a.chip)
    inp = pyBigWig.open(a.input) if a.input else None
    for bw, path in ((chip, a.chip), (inp, a.input)):
        if bw is not None and a.chrom not in bw.chroms():
            raise SystemExit(f"{path} has no {a.chrom}")

    def background(chrom, mid):
        if inp is None:
            return 0.0
        return mean_signal(inp, chrom, mid - a.input_window // 2, mid + a.input_window // 2)

    peaks = {}
    with open(a.peaks) as fh:
        for line in fh:
            r = line.rstrip("\n").split("\t")
            if line.startswith(("#", "track", "browser")) or len(r) < 4 or r[0] != a.chrom:
                continue
            s, e = int(r[1]), int(r[2])
            off = int(r[9]) if len(r) > 9 and int(r[9]) >= 0 else (e - s) // 2
            summit = s + off
            S = mean_signal(chip, a.chrom, summit - a.window // 2, summit + a.window // 2)
            I = background(a.chrom, summit)
            peaks[r[3]] = (S, I, max(S - I, 0.0))
    if not peaks:
        raise SystemExit(f"no peaks on {a.chrom} in {a.peaks}")
    E = np.array([p[2] for p in peaks.values()])
    ref = float(E.max()) if a.ref == "max" else float(np.percentile(E, float(a.ref.lstrip("p"))))
    if ref <= 0:
        raise SystemExit("reference enrichment is 0: ChIP is not above input at any peak")

    rows, missing = [], 0
    with open(a.sites) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            r = line.rstrip("\n").split("\t")
            name = r[3]
            if name in peaks:
                S, I, e = peaks[name]
            else:                              # site without its peak: measure it directly
                missing += 1
                mid = (int(r[1]) + int(r[2])) // 2
                S = mean_signal(chip, r[0], mid - a.window // 2, mid + a.window // 2)
                I = background(r[0], mid)
                e = max(S - I, 0.0)
            rows.append((name, r[0], r[1], r[2], r[5] if len(r) > 5 else ".", S, I, e, min(1.0, e / ref)))

    occ = np.array([r[-1] for r in rows]) if rows else np.zeros(0)
    with open(a.output, "w") as out:
        out.write(f"# chip={a.chip} input={a.input or 'none (ChIP alone)'} chrom={a.chrom} window={a.window} "
                  f"input_window={a.input_window}\n")
        out.write(f"# chromosome peaks={len(peaks)} ref({a.ref})={ref:.4f} "
                  f"max_enrichment={E.max():.4f} median_enrichment={np.median(E):.4f}\n")
        out.write("name\tchrom\tstart\tend\tstrand\tchip\tinput\tenrichment\toccupancy\n")
        for r in rows:
            out.write("\t".join(map(str, r[:5])) + "\t" + "\t".join(f"{v:.5g}" for v in r[5:]) + "\n")
    print(f"{len(rows)} sites ({missing} measured without a matching peak); ref({a.ref})={ref:.3f}; "
          f"occupancy median {np.median(occ) if len(occ) else float('nan'):.3f}, "
          f"max {occ.max() if len(occ) else float('nan'):.3f}")


if __name__ == "__main__":
    main()
