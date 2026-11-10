#!/usr/bin/env python
"""
Ebf3 Stage 3: enhancer tiers and the Ebf3 promoter as 10-column narrowPeak files
for ChromSimPipe's `sticky:` option, plus a per-element evidence table.

Evidence classes (kept separate in the table):
  annotation  VISTA element (UCSC mm10 vistaEnhancers.bb) / GENCODE vM25 TSS
  peak        overlap with a GSE96107 IDR-0.05 narrowPeak
  signal      mean bigWig signal over the element / mean Input over the element
              (GSE96107 bigWigs; a descriptive ratio, NOT a statistical test)

Tiers (separate files, never mixed):
  tier1  hs737 + (Ebf3 promoter, separate file)
  tier2  VISTA elements in the window with CN H3K27ac/Input > 2
  tier3  all IDR-0.05 H3K27ac peaks of the cell type in the window, minus peaks
         within --promoter-flank of a GENCODE TSS of Mgmt, Glrx3, 9430038I01Rik
         or Ebf3. ES has no H3K27ac peaks in GSE96107 -> no ES tier3.
Tier1/2 element SETS are defined once (from annotation + CN signal) and written
for every cell type; only column 7 (signalValue) changes per cell type.

The "summit" of a VISTA element (which has none) is its midpoint (column 10 =
midpoint - start).

Run with the ebf3_analysis env (pyBigWig, bigBedToBed), from the repo root:
  python -I scripts/ebf3/build_regulatory_sites.py --window chr7:135800000-137700000
"""
import argparse, gzip, json, os, subprocess, sys
import numpy as np
import pyBigWig

RAW = "data/raw"
GEO = f"{RAW}/gse96107"
EBF3_TSS = 137_314_445          # GENCODE vM25 / RefSeq NM_001113415.1 txEnd (minus strand)
EXCLUDE_GENES = ("Mgmt", "Glrx3", "9430038I01Rik", "Ebf3")
TIER2_CUTOFF = 2.0


def region(s):
    c, r = s.replace(",", "").split(":")
    a, b = r.split("-")
    return c, int(a), int(b)


def vista(chrom, start, end):
    out = subprocess.run(["bigBedToBed", f"-chrom={chrom}", f"-start={start}", f"-end={end}",
                          f"{RAW}/ucsc/vistaEnhancers.bb", "stdout"],
                         check=True, capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        r = line.split("\t")
        rows.append(dict(id=r[3], chrom=r[0], start=int(r[1]), end=int(r[2]), raw=r))
    return rows


def peaks(name, chrom, start, end):
    path = f"{GEO}/GSE96107_{name}.IDR0.05.filt.narrowPeak.gz"
    if not os.path.exists(path):
        return None
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            r = line.rstrip("\n").split("\t")
            if r[0] == chrom and int(r[2]) > start and int(r[1]) < end:
                out.append(r)
    return sorted(out, key=lambda r: int(r[1]))


def tss_of(genes, chrom):
    tss = {g: set() for g in genes}
    with gzip.open(f"{RAW}/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            r = line.split("\t")
            if r[0] != chrom or r[2] != "transcript":
                continue
            g = r[8].split('gene_name "')[1].split('"')[0]
            if g in tss:
                tss[g].add(int(r[3]) if r[6] == "+" else int(r[4]))
    return tss


def all_tss(chrom, start, end):
    out = []
    with gzip.open(f"{RAW}/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            r = line.split("\t")
            if r[0] != chrom or r[2] != "transcript":
                continue
            t = int(r[3]) if r[6] == "+" else int(r[4])
            if start <= t < end:
                out.append((t, r[8].split('gene_name "')[1].split('"')[0]))
    return out


class Signal:
    def __init__(self):
        self.bw = {}

    def mean(self, name, chrom, s, e):
        path = f"{GEO}/GSE96107_{name}.bw"
        if not os.path.exists(path):
            return float("nan")
        if name not in self.bw:
            self.bw[name] = pyBigWig.open(path)
        v = self.bw[name].stats(chrom, s, e, type="mean", exact=True)[0]
        return float(v) if v is not None else 0.0

    def ratio(self, mark, cell, chrom, s, e):
        inp = self.mean(f"{cell}_Input", chrom, s, e)
        m = self.mean(f"{cell}_{mark}", chrom, s, e)
        return m / inp if inp > 0 else float("nan")


def np_row(chrom, s, e, name, signal, summit_bp, score=1000):
    sig = 0.0 if not np.isfinite(signal) else signal
    return f"{chrom}\t{s}\t{e}\t{name}\t{score}\t.\t{sig:.5f}\t-1\t-1\t{summit_bp - s}\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="chr7:135800000-137700000")
    ap.add_argument("--annot-window", default="chr7:135400000-138100000",
                    help="wider window for the evidence table (reports mm1590 etc.)")
    ap.add_argument("--promoter-flank", type=int, default=1000)
    ap.add_argument("--resolution", type=int, default=1000)
    ap.add_argument("--outdir", default="data/peaks")
    ap.add_argument("--table", default="results/ebf3/regulatory_evidence.tsv")
    a = ap.parse_args()
    chrom, ws, we = region(a.window)
    _, aws, awe = region(a.annot_window)
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(os.path.dirname(a.table), exist_ok=True)
    sig = Signal()
    mon = lambda bp: (bp - ws) // a.resolution

    # ---------------- evidence table over the wide window ----------------
    pk = {f"{c}_{m}": peaks(f"{c}_{m}", chrom, aws, awe)
          for c in ("ES", "NPC", "CN") for m in ("H3K27ac", "H3K4me3", "CTCF")}
    rows = []
    for v in vista(chrom, aws, awe):
        mid = (v["start"] + v["end"]) // 2
        rec = dict(element=v["id"], chrom=chrom, start=v["start"], end=v["end"], summit=mid,
                   in_sim_window=ws <= mid < we, monomer=mon(mid) if ws <= mid < we else "",
                   evidence="annotation(VISTA)")
        for c in ("NPC", "CN"):
            for m in ("H3K27ac", "H3K4me1", "H3K4me3"):
                rec[f"{c}_{m}_over_input"] = round(sig.ratio(m, c, chrom, v["start"], v["end"]), 3)
        for k, p in pk.items():
            rec[f"{k}_IDRpeak"] = ("NA" if p is None else
                                   int(any(int(r[2]) > v["start"] and int(r[1]) < v["end"] for r in p)))
        rec["accessibility_proxy"] = "none available in GSE96107 (no ATAC/DNase)"
        rows.append(rec)
    rows.sort(key=lambda r: r["start"])

    # ---------------- tiers ----------------
    in_win = [r for r in rows if r["in_sim_window"]]
    hs737 = [r for r in in_win if r["element"] == "hs737"]
    if len(hs737) != 1:
        sys.exit(f"hs737 not found exactly once in {a.window}: {[r['element'] for r in in_win]}")
    tier1 = hs737
    tier2 = [r for r in in_win if r["CN_H3K27ac_over_input"] > TIER2_CUTOFF]
    for r in rows:
        r["tier"] = ",".join(t for t, members in (("1", tier1), ("2", tier2)) if r in members) or "-"

    tss = tss_of(EXCLUDE_GENES, chrom)
    excl = sorted(t for g in EXCLUDE_GENES for t in tss[g])
    other_tss = [(t, g) for t, g in all_tss(chrom, ws, we) if g not in EXCLUDE_GENES]

    def near(s, e, points, flank):
        return [p for p in points if s - flank <= (p[0] if isinstance(p, tuple) else p) < e + flank]

    summary = {"window": a.window, "promoter_flank_bp": a.promoter_flank,
               "excluded_TSS": {g: sorted(tss[g]) for g in EXCLUDE_GENES}, "files": {}}
    for cell in ("ES", "NPC", "CN"):
        for tier, members in (("tier1", tier1), ("tier2", tier2)):
            path = f"{a.outdir}/{cell}_enhancers_{tier}.narrowPeak"
            with open(path, "w") as f:
                for r in members:
                    s = sig.ratio("H3K27ac", cell, chrom, r["start"], r["end"]) if cell != "ES" else float("nan")
                    f.write(np_row(chrom, r["start"], r["end"], r["element"], s, r["summit"]))
            summary["files"][path] = [r["element"] for r in members]
        ac = pk[f"{cell}_H3K27ac"]
        if ac is None:
            summary["files"][f"{a.outdir}/{cell}_enhancers_tier3.narrowPeak"] = \
                f"NOT WRITTEN: GSE96107 has no {cell} H3K27ac IDR peaks"
            continue
        kept, dropped, on_other = [], [], []
        for r in ac:
            s, e = int(r[1]), int(r[2])
            summit = s + int(r[9])
            if not ws <= summit < we:
                continue
            (dropped if near(s, e, excl, a.promoter_flank) else kept).append(r)
        path = f"{a.outdir}/{cell}_enhancers_tier3.narrowPeak"
        with open(path, "w") as f:
            for i, r in enumerate(kept, 1):
                s = int(r[1])
                hit = near(s, int(r[2]), other_tss, a.promoter_flank)
                if hit:
                    on_other.append(f"{cell}_H3K27ac_{i:03d}@{s + int(r[9])}:" + "/".join(sorted({g for _, g in hit})))
                f.write("\t".join([chrom, r[1], r[2], f"{cell}_H3K27ac_{i:03d}", r[4], ".",
                                   r[6], r[7], r[8], r[9]]) + "\n")
        summary["files"][path] = dict(n=len(kept), dropped_at_excluded_promoters=len(dropped),
                                      kept_but_near_other_annotated_TSS=on_other)

    # ---------------- promoter ----------------
    prom = f"{a.outdir}/Ebf3_promoter.narrowPeak"
    ps, pe = EBF3_TSS - 500, EBF3_TSS + 500
    with open(prom, "w") as f:
        f.write(np_row(chrom, ps, pe, "Ebf3_TSS", sig.ratio("H3K27ac", "CN", chrom, ps, pe), EBF3_TSS))
    summary["files"][prom] = dict(tss=EBF3_TSS, summit_monomer=mon(EBF3_TSS),
                                  peak_monomers=[mon(ps), mon(pe - 1)])
    rows.append(dict(element="Ebf3_promoter", chrom=chrom, start=ps, end=pe, summit=EBF3_TSS,
                     in_sim_window=True, monomer=mon(EBF3_TSS),
                     evidence="annotation(GENCODE vM25 TSS)",
                     **{f"{c}_{m}_over_input": round(sig.ratio(m, c, chrom, ps, pe), 3)
                        for c in ("NPC", "CN") for m in ("H3K27ac", "H3K4me1", "H3K4me3")},
                     **{f"{k}_IDRpeak": ("NA" if p is None else
                                         int(any(int(r[2]) > ps and int(r[1]) < pe for r in p)))
                        for k, p in pk.items()},
                     accessibility_proxy="none available in GSE96107 (no ATAC/DNase)", tier="promoter"))

    # E and P summits sharing a bin
    e_mon = {r["element"]: mon(r["summit"]) for r in tier1 + tier2}
    summary["E_P_shared_bin"] = [k for k, m in e_mon.items() if m == mon(EBF3_TSS)]

    cols = list(rows[0].keys())
    with open(a.table, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(a.table.replace(".tsv", "_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
