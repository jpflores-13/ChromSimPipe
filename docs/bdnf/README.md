# Bdnf locus: design and caveats

Config `config/ChromSimConfig.bdnf.yaml`, samplesheet `ChromSimSamplesheet.bdnf.txt`,
output `output_bdnf/`. Same five models (A, A0, B, C, D), polymer, cohesin settings and
sampling as the Ebf3 run; separate config so the Ebf3 inputs never change.

## Locus

Bdnf, GENCODE vM25: chr2:109,674,700–109,727,007 (+), 9 alternative TSSs (109,674,700;
109,675,897; 109,677,032; 109,692,436; 109,693,027; 109,693,563; 109,694,409;
109,694,469; 109,709,383). Window **chr2:108,800,000–110,300,000** (1.5 Mb, 1,500
monomers, 34 tiled copies). Edges chosen in CTCF-free gaps (ES/CN IDR peaks): left
between the 108.67–108.71 Mb cluster and 108.88 Mb, right between 110.23 and 110.36 Mb.
Insulation edge check (same ±20 kb / chromosome-p75 rule as Ebf3) runs in the gate job
before any GPU job: `results/bdnf/window_edge_check.json`.

## Caveats (read before interpreting)

1. **Bdnf is silent in every Bonev 2017 state**: 0.03–0.62 FPKM in all 10 GSE96107
   RNA samples (ES, in vitro NPC/CN, ncx_NPC/CN). No "on" state exists in this dataset.
2. **No validated distal enhancer.** No VISTA element within 1.4 Mb. The validated
   Bdnf enhancer is intronic and acts on the exon-I transcripts (Tuvikene et al. 2021,
   eLife, doi:10.7554/eLife.65161); its exact mm10 coordinates are not yet taken from
   the paper. Two candidate summits 1–4 kb downstream of the exon-I TSS cluster
   (109,678,118; 109,678,648) may correspond to it; unverified.
3. **The sticky enhancer set is permissive and mostly a proxy**: 26 summits on 16
   monomers = union of CN H3K27ac IDR summits (3) and ENCODE ENCSR810HQR E14.5
   forebrain ATAC IDR summits (23, PROXY), minus summits within ±1 kb of any
   protein-coding TSS. Promoters = all 9 Bdnf TSSs (7 monomers). E–E = P–P = 0,
   E–P = 3 kT (verified with `scripts/ebf3/check_sticky_types.py`). Many candidates
   lie near Lin7c (109.90–109.96 Mb). C and D are an exploratory upper bound, not a
   test of a validated E–P pair. Condition names are kept identical to Ebf3
   (`C_sticky_t1`, `D_combined_t1`) so figures share code and colours; "t1" does not
   mean tier 1 here.
4. NPC / CN Hi-C is neocortex (ncx_); CTCF / H3K27ac are in vitro (same mismatch as Ebf3).

Scripts: `scripts/ebf3/build_bdnf_sites.py`, `scripts/ebf3/window_edge_check.py`,
`scripts/ebf3/bdnf_gate.sh`.
