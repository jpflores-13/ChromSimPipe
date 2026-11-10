# Ebf3 Stage 2: CTCF barriers

## Inputs

GSE96107 (Bonev et al. 2017) IDR-0.05 CTCF peaks, in vitro ES / NPC / CN:
`GSE96107_{ES,NPC,CN}_CTCF.IDR0.05.filt.narrowPeak.gz` (in
`data/raw/gse96107/`). Genome-wide: ES 121,406, NPC 34,367, CN 28,146.

Working copies `data/peaks/{ES,NPC,CN}_CTCF.narrowPeak` are the chr7 rows,
sorted, with column 4 replaced by a unique name (`CN_CTCF_chr7_00001`, ...).
Everything else is unchanged. **Why:** every GSE96107 peak is named `.`, and
`configs/parameters.get_ctcf_signals` matches oriented sites to their peak
signal by name. With the original names all sites receive one signal, and
`ctcf_capture: {min, max}` silently becomes a uniform `(min + max) / 2`.

Orientation: `scripts/extract_ctcf_sites.py` with the config's `ctcf_motif`
(`MA0139.1`, FIMO p < 1e-4, ±250 bp around the summit, `unoriented: drop`),
genome `data/mm10.primary.fa`. The workflow's `orientCTCF` orients only the
`ctcf_type`s used by a condition (CN), so ES and NPC were run with the same
arguments by hand. Output: `output_ebf3/loci/ebf3_hs737/ctcf/ctcf_oriented_mm10_{ES,NPC,CN}.bed`
(logs alongside in `output_ebf3/logs/ebf3_hs737/`). Peaks are fixed-width with
the summit at the centre, so the BED midpoint that `load_ctcf_from_bed` uses
equals the summit.

`scripts/validate_ctcf.py` (CN, the type used by every condition): 27 passed,
0 warnings, 0 failed.

## Counts in chr7:135,800,000–137,700,000

| Cell type | Peaks in window | Oriented by FIMO | Dropped (no hit) | + / − | Convergent (+ … −) pairs, all | Adjacent convergent pairs | Convergent pairs enclosing both hs737 and the Ebf3 TSS |
|---|---|---|---|---|---|---|---|
| ES | 82 | 67 | 15 | 37 / 30 | 610 | 16 | 72 |
| NPC | 24 | 20 | 4 | 15 / 5 | 43 | 5 | 6 |
| CN | 13 | 12 | 1 | 10 / 2 | 16 | 2 | 6 |

Peak counts match the brief's expected 82 / 24 / 13. Dropped peaks (summits):
ES 135,853,664; 135,968,828; 136,021,089; 136,055,119; 136,121,596;
136,136,459; 136,209,472; 136,333,004; 136,537,243; 136,879,113; 136,938,112;
137,122,439; 137,308,523; 137,406,300; 137,484,344. NPC 135,853,663;
136,055,068; 136,121,589; 137,471,937. CN 135,968,801.

CN sites (summit kb, strand, signalValue): 135,837.2 + (184); 135,889.5 + (37);
136,010.4 + (62); 136,152.7 + (266); 136,294.6 + (71); 136,489.1 + (320);
137,246.7 + (191); **137,321.7 − (207)**; 137,468.4 + (212); 137,478.4 + (476);
137,550.2 − (132); 137,603.3 + (244). Machine-readable: `results/ebf3/ctcf_barriers.json`.

All eleven candidate anchors in the brief (135,837 / 136,010 / 136,152 /
136,294 / 136,489 and 137,246 / 137,321 / 137,468 / 137,478 / 137,550 /
137,603 kb) are present and oriented in CN; all are also present in ES.

## Observations relevant to the E–P question (observation, not interpretation)

- In CN the only − site right of the gene body that could stop a rightward leg
  near the promoter is at 137,321.7 kb, **7 kb downstream of the Ebf3 TSS**
  (monomer 1521 vs promoter monomer 1514). Left of hs737 (monomer 283) are +
  sites at 136,010 (m210) and, to its right, 136,152 (m352). A convergent pair
  136,010(+) … 137,321.7(−) therefore nearly brackets the E–P pair. Loop
  extrusion with CTCF stalling alone (model B) may thus raise contact
  frequency close to the hs737 × promoter pixel. The E–P endpoint window must
  be small enough to separate the E–P pixel from the CTCF corner at
  (136,010 kb, 137,321.7 kb), which is 73 kb (left) and 7 kb (right) away.
- CN has only 2 − sites in the window. There is no − site at the left of the
  window, so leftward-moving legs left of 135,837 kb run to the window edge.

## Two caveats (flagged, not modelled)

1. **CN has 4.3× fewer CTCF peaks than ES genome-wide (28,146 vs 121,406).**
   The most likely explanation is peak-call sensitivity (depth, IP efficiency,
   IDR on fewer replicates), not loss of CTCF binding. The ES → CN barrier
   difference in the window (67 → 12 oriented sites) must not be interpreted
   as biology without renormalisation, e.g. re-calling peaks from the raw
   reads at matched depth and threshold, or quantifying CTCF bigWig signal at
   the union of sites. The `ctcf_type: ES vs CN` sensitivity condition tests
   the model's sensitivity to the barrier set, not a biological ES → CN change.
2. **`ctcf_capture: {min, max}` ranks sites by signalValue within the locus**
   (`sim_features.capture_by_signal`). With 12 CN sites the rank spread is
   coarse, and each site's capture is fixed by its rank among 12. `min` and
   `max` are free parameters on a sensitivity axis, not measurements.

## Figure

`figures/ebf3/03_ctcf_orientation.png` (+ `.pdf`): strand-coded arrows per cell
type, adjacent convergent pairs bracketed (all-pairs would be 610 brackets in
ES), hs737 and the Ebf3 TSS marked and linked by an arc, the conserved anchors
as keyed triangles. Script: `scripts/ebf3/plot_ctcf_orientation.py`.
