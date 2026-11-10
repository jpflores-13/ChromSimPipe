# Ebf3 Stage 4: configuration rationale

Files: `config/ChromSimConfig.ebf3.yaml` (production),
`config/ChromSimConfig.ebf3.test.yaml` (complete copy: 1 replicate, 1 shard,
2,000 frames, `outputDir: output_test_ebf3`; nothing else differs; checked
with `diff`), `ChromSimSamplesheet.ebf3.txt`. Both configs parse with
`configs/parameters.py`. Every key used exists in `configs/parameters.py` at
`12d5e80`.

## Global settings

| Setting | Value | Why |
|---|---|---|
| `samplesheet` | `ChromSimSamplesheet.ebf3.txt` | separate from the STRS samplesheet; the shipped one is untouched |
| `fileNamesFrom` | `['Condition']` | sample names ES / NPC / CN |
| `outputDir` | `output_ebf3` | own tree; `output_*` is git-ignored |
| `genome`, `genome_assembly` | `data/mm10.primary.fa`, `mm10` | UCSC mm10 primary chromosomes, md5-checked (00_environment_audit.md) |
| `ctcf_motif` | MA0139.1, p < 1e-4, ±250 bp, drop | brief; MA0139.1 is the vertebrate CTCF core motif, conserved in mouse |
| `resolution` | 1000 | 1 monomer = 1 kb = the finest Hi-C bin |
| `n_replicates`, `n_shards` | 3, 4 | brief (production) |
| `hpcProfile` | `profiles/longleaf` | verified against `sacctmgr` / `sinfo` |
| `analysis.hic_balance` | `weight` | the mcools are ICE-balanced by `cooler zoomify --balance` (cis-only); there is no KR column. `auto` would choose the same; explicit is clearer |
| `tiling` | padding 270, chain 70,000 | shipped values; 1,900 + 540 = 2,440 per copy → **28 copies**, 68,320-monomer chain |
| `sim_run.frames_per_replicate` | 10,000 | shipped value. 2,500 frames/shard × `saveevery` 50 = 125,000 blocks, a multiple of `blocks_per_updater` 1000 (required). Test: 2,000 × 50 = 100,000 ✓ |
| `polysim:` | shipped bootcamp values, incl. `initskip: 20` | unchanged; identical across all conditions so A/A0/B/C/D are matched |

## Locus

`ebf3_hs737`, chr7:135,800,000–137,700,000 (1,900 monomers), as specified:
283 kb left of hs737, 386 kb right of the Ebf3 TSS (both > 250 kb), encloses
the conserved CTCF clusters, edges on the 1 kb grid. `anchor_bp: 137314000` =
the promoter bin (monomer 1514); it only affects the optional MSD probes (MSD
is off). The window-edge insulation check is in 01_hic_qc.md; if it fails the
window becomes 135,467,000–137,930,000 and every downstream file is regenerated.

## Samplesheet

| Column | Notes |
|---|---|
| `CTCF_Peaks` | `data/peaks/{ES,NPC,CN}_CTCF.narrowPeak`: chr7 rows of the GSE96107 IDR files, renamed uniquely (02_ctcf_barriers.md) |
| `Hi-C` | ES → `ES_chr7.mcool`; **NPC → `ncx_NPC_chr7.mcool`, CN → `ncx_CN_chr7.mcool`** (neocortex Hi-C paired with in vitro ChIP: mismatch, see 01_hic_qc.md) |
| `Enhancer_Peaks` | tier 1 (what the brief's conditions name) |
| `Enhancer_Peaks_T2`, `Enhancer_Peaks_T3` | **added**: `sticky.peaks` names a samplesheet column, so each tier needs its own column. ES tier 3 = `NA` (no ES H3K27ac in GSE96107) |
| `Promoter_Peaks` | `data/peaks/Ebf3_promoter.narrowPeak` (all rows) |

## Conditions and the comparisons they support

| Condition | params_name | Cohesin | CTCF | Polymer force | Sticky | Compared with |
|---|---|---|---|---|---|---|
| A_polymer | polymer | none | (ignored) | default repulsion | none | B |
| A0_polymer_attr | polymer_attr | none | (ignored) | heteropolymer_SSW, `attraction_energy: 0` | none | C, D |
| B_extrusion | extrusion | lifetime 75,000, sep 240 | CN, capture 0.5 | default | none | A |
| C_sticky_t1 | sticky_t1 | none | (ignored) | heteropolymer_SSW | E–P 3 kT | A0 |
| D_combined_t1 | combined_t1 | lifetime 75,000, sep 240 | CN, capture 0.5 | heteropolymer_SSW | E–P 3 kT | A0 (and B, with the force caveat) |

**Why A0 is mandatory.** `run_simulation_shard.py` passes `monomer_types` /
`interaction_matrix` to polysim when a condition has `sticky`, and when it sets
`attraction_energy` (even 0) it passes one type with a zero matrix. Either
switches polysim to its heteropolymer force, whose repulsion has a different
shape from the default. C − A would mix mechanism with force model; C − A0 and
D − A0 do not. **D − B also crosses force models**; read it only alongside
D − A0 and B − A.

Cohesin values are the polysim bootcamp defaults (processivity ≈ 75,000 ×
0.0025 ≈ 190 kb per leg; one cohesin per 240 kb). They are starting points to
be calibrated on ES P(s) and boundaries (Stage 6), not fitted values.

### The 5 vs 6 condition count

The brief's condition list has 5 entries (A, A0, B, C_t1, D_t1) but its
production budget says 6 conditions (72 jobs). The config has the 5 listed (=
60 GPU jobs). The sixth is not specified and was not invented.

## Sticky sites: verification of the interaction matrix

`energy: 0.0` per factor plus `sticky_cross` gives attraction only between
enhancers and promoters. A single factor holding both would make every enhancer
attract every other enhancer (and every promoter). Output of
`scripts/ebf3/check_sticky_types.py` (the same calls `run_simulation_shard.py`
makes), for `config/ChromSimConfig.ebf3.yaml`:

```
== C_sticky_t1  (chain 68320 monomers, 28 copies)
  Enhancer: local monomers [283] -> bp [136083000]
  Promoter: local monomers [1514] -> bp [137314000]
  E and P summits sharing a bin: none
  labels: ['none', 'Enhancer', 'Promoter']
  matrix (kT):
[[0. 0. 0.]
 [0. 0. 3.]
 [0. 3. 0.]]
  type counts on the chain: {0: 68264, 1: 28, 2: 28}
  types in copy 0 at those monomers: {283: 1, 1514: 2}
== D_combined_t1  (chain 68320 monomers, 28 copies)
  (identical)
```

E–E = 0, P–P = 0, E–P = 3.0 kT, as intended. No monomer holds both E and P, so
no combined self-attracting type exists for tier 1 (also none for tier 2).
`_parse_sticky` imposes no positivity constraint, so `energy: 0.0` is accepted.
Tier 3: the E site at monomer 1509 is 5 monomers from P (1514), but they are
different monomers, so no combined type.

### Modelling caveat found during verification

Types are global on the tiled chain, so **hs737 in one copy attracts the Ebf3
promoter of every other copy** (28 of each on one confined chain, copies
separated by 540 padding monomers). Inter-copy E–P contacts never appear in
the averaged per-copy contact map, but they can occupy sticky sites and so
*lower* the intra-copy E–P frequency. At 3 kT with one E and one P per copy
this is likely small but unmeasured. Check: count inter-copy vs intra-copy E–P
contacts in the C and D trajectories (planned in Stage 6).

## Model E (sensitivity), one key changed per condition

Not in the production config until the first batch's wall time is measured
(Stage 5). Base for cohesin axes: B; base for sticky axes: D (and C for energy).

| Axis | Values | Base condition |
|---|---|---|
| `ctcf_capture` | 0.25, 0.5, 0.9, {min 0.2, max 0.9} | B |
| `lifetime` | 37,500, 75,000, 150,000 | B |
| `separation` | 120, 240, 480 | B |
| sticky E–P energy | 1, 2, 3, 4 kT | C and D |
| tier | 1, 2, 3 (columns `Enhancer_Peaks`, `_T2`, `_T3`) | D |
| `span` | summit vs peak | D |
| `ctcf_type` | ES vs CN, cohesin fixed | B (barrier-set sensitivity only; see 02_ctcf_barriers.md caveat 1) |

`params_name` must be unique within a `ctcf_type`; every sweep entry gets its
own.
