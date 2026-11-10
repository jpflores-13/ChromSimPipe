#!/usr/bin/env python
"""
Bdnf locus: promoter and candidate-enhancer sticky sites (10-column narrowPeak).

Bdnf has NO validated distal enhancer (no VISTA element within 1.4 Mb; the
validated enhancer is intronic, Tuvikene et al. 2021 eLife,
doi:10.7554/eLife.65161) and is silent in every Bonev 2017 state (<= 0.62 FPKM).
The enhancer set here is therefore PERMISSIVE and exploratory (like Ebf3 tier 3):

  promoters  every GENCODE vM25 Bdnf transcript TSS (9 alternative promoters), TSS +- 500 bp
  enhancers  union of CN H3K27ac IDR-0.05 summits (GSE96107, in vitro CN) and
             ATAC IDR summits (ENCODE ENCSR810HQR, E14.5 forebrain; PROXY) in the window,
             minus summits within +-1 kb of ANY protein-coding TSS in the window

Outputs: data/peaks/bdnf/{CN,NPC,ES}_enhancers.narrowPeak (same element set; ES/NPC rows
identical, the set is defined once), data/peaks/bdnf/Bdnf_promoters.narrowPeak,
results/bdnf/sites_summary.json
Run from the repo root with the ebf3_analysis env.
"""
import gzip, json, os

CHROM, WS, WE = "chr2", 108_800_000, 110_300_000
FLANK = 1000
OUT = "data/peaks/bdnf"
os.makedirs(OUT, exist_ok=True)
os.makedirs("results/bdnf", exist_ok=True)


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


def summits(path, label):
    out = {}
    with gzip.open(path, "rt") as f:
        for l in f:
            r = l.rstrip("\n").split("\t")
            if r[0] != CHROM:
                continue
            s = int(r[1]) + int(r[9])
            if WS <= s < WE and (s not in out or float(r[6]) > out[s][1]):
                out[s] = (label, float(r[6]))
    return out


cand = summits("data/raw/gse96107/GSE96107_CN_H3K27ac.IDR0.05.filt.narrowPeak.gz", "CN_H3K27ac")
for s, v in summits("data/raw/encode_atac/ENCFF776MOH.bed.gz", "ATAC_forebrainE14.5").items():
    cand.setdefault(s, v)
kept = {s: v for s, v in sorted(cand.items()) if not any(abs(s - t) <= FLANK for t in pc_tss)}

for cell in ("CN", "NPC", "ES"):
    with open(f"{OUT}/{cell}_enhancers.narrowPeak", "w") as f:
        for i, (s, (lab, sig)) in enumerate(kept.items(), 1):
            f.write(f"{CHROM}\t{s - 250}\t{s + 250}\tBdnfE_{i:03d}_{lab}\t1000\t.\t{sig:.5f}\t-1\t-1\t250\n")
with open(f"{OUT}/Bdnf_promoters.narrowPeak", "w") as f:
    for i, t in enumerate(bdnf_tss, 1):
        f.write(f"{CHROM}\t{t - 500}\t{t + 500}\tBdnf_TSS{i}\t1000\t.\t1.0\t-1\t-1\t500\n")

mon = lambda bp: (bp - WS) // 1000
e_mon, p_mon = sorted({mon(s) for s in kept}), sorted({mon(t) for t in bdnf_tss})
summary = dict(window=f"{CHROM}:{WS}-{WE}", bdnf_tss=bdnf_tss, promoter_monomers=p_mon,
               candidates=len(cand), excluded_at_protein_coding_TSS=len(cand) - len(kept),
               enhancers=len(kept), enhancer_monomers=len(e_mon),
               by_source={k: sum(1 for v in kept.values() if v[0] == k) for k in ("CN_H3K27ac", "ATAC_forebrainE14.5")},
               shared_E_P_monomers=sorted(set(e_mon) & set(p_mon)),
               E_within_5kb_of_a_promoter=[s for s in kept if any(abs(s - t) <= 5000 for t in bdnf_tss)])
json.dump(summary, open("results/bdnf/sites_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
