#!/usr/bin/env python
"""
Bdnf extra-wide window: promoter and candidate-enhancer sticky sites, with the two
enhancer evidence types kept as SEPARATE files / sticky factors (build_bdnf_sites.py merged
them into one deduplicated union). The config uses both as enhancers in the same simulation.

  H3K27ac  CN H3K27ac IDR-0.05 summits (GSE96107, in vitro CN)
  ATAC     ENCODE ENCSR810HQR E14.5 forebrain ATAC IDR summits (PROXY, not a cell-type match)
  promoter every GENCODE vM25 Bdnf transcript TSS +- 500 bp (same 9 TSSs)

Each set drops summits within +-1 kb of ANY protein-coding TSS in the window. A summit
present in both sets stays in both (monomers holding both factors get their own type).

Outputs (10-column narrowPeak, same rows for ES/NPC/CN, like the original set):
  data/peaks/bdnf_xwide/{CN,NPC,ES}_enhancers_{H3K27ac,ATAC}.narrowPeak
  data/peaks/bdnf_xwide/Bdnf_promoters.narrowPeak
  results/bdnf_xwide/sites_summary.json
Run from the repo root with the ebf3_analysis env:
  .conda_envs/ebf3_analysis/bin/python -I scripts/ebf3/build_bdnf_xwide_sites.py
"""
import gzip, json, os

CHROM, WS, WE = "chr2", 106_580_000, 112_100_000
FLANK = 1000
OUT, RES = "data/peaks/bdnf_xwide", "results/bdnf_xwide"
os.makedirs(OUT, exist_ok=True)
os.makedirs(RES, exist_ok=True)

SOURCES = {
    "H3K27ac": "data/raw/gse96107/GSE96107_CN_H3K27ac.IDR0.05.filt.narrowPeak.gz",
    "ATAC":    "data/raw/encode_atac/ENCFF776MOH.bed.gz",
}


def gtf():
    with gzip.open("data/raw/gencode/gencode.vM25.annotation.gtf.gz", "rt") as f:
        for l in f:
            if not l.startswith("#"):
                yield l.rstrip("\n").split("\t")


bdnf_tss, pc_tss = set(), set()
for r in gtf():
    if r[0] != CHROM or r[2] != "transcript":
        continue
    t = int(r[3]) if r[6] == "+" else int(r[4])
    if 'gene_name "Bdnf"' in r[8]:
        bdnf_tss.add(t)
    if 'gene_type "protein_coding"' in r[8] and WS <= t < WE:
        pc_tss.add(t)
bdnf_tss = sorted(bdnf_tss)


def summits(path):
    out = {}
    with gzip.open(path, "rt") as f:
        for l in f:
            r = l.rstrip("\n").split("\t")
            if r[0] != CHROM:
                continue
            s = int(r[1]) + int(r[9])
            if WS <= s < WE and (s not in out or float(r[6]) > out[s]):
                out[s] = float(r[6])
    return out


mon = lambda bp: (bp - WS) // 1000
summary = dict(window=f"{CHROM}:{WS}-{WE}", bdnf_tss=bdnf_tss,
               promoter_monomers=sorted({mon(t) for t in bdnf_tss}), sets={})
for name, path in SOURCES.items():
    cand = summits(path)
    kept = {s: v for s, v in sorted(cand.items()) if not any(abs(s - t) <= FLANK for t in pc_tss)}
    for cell in ("CN", "NPC", "ES"):
        with open(f"{OUT}/{cell}_enhancers_{name}.narrowPeak", "w") as f:
            for i, (s, sig) in enumerate(kept.items(), 1):
                f.write(f"{CHROM}\t{s - 250}\t{s + 250}\tBdnfE_{name}_{i:03d}\t1000\t.\t{sig:.5f}\t-1\t-1\t250\n")
    summary["sets"][name] = dict(
        candidates=len(cand), excluded_at_protein_coding_TSS=len(cand) - len(kept),
        enhancers=len(kept), enhancer_monomers=len({mon(s) for s in kept}),
        shared_E_P_monomers=sorted({mon(s) for s in kept} & set(summary["promoter_monomers"])),
        n_left_of_first_TSS=sum(1 for s in kept if s < bdnf_tss[0]),
        n_right_of_last_TSS=sum(1 for s in kept if s > bdnf_tss[-1]))

with open(f"{OUT}/Bdnf_promoters.narrowPeak", "w") as f:
    for i, t in enumerate(bdnf_tss, 1):
        f.write(f"{CHROM}\t{t - 500}\t{t + 500}\tBdnf_TSS{i}\t1000\t.\t1.0\t-1\t-1\t500\n")
json.dump(summary, open(f"{RES}/sites_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
