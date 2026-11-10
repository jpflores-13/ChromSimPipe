#!/usr/bin/env python
"""
Build docs/ebf3/Ebf3_report.md and render docs/ebf3/Ebf3_report.pdf (pandoc +
xelatex). Sections whose inputs do not exist yet (simulation figures, endpoint
table) are replaced by a "pending" note, so this can run at any stage.
Run from the repo root with the ebf3_analysis env:
  python -I scripts/ebf3/build_report.py
"""
import datetime, json, os, subprocess
import pandas as pd

OUT_MD, OUT_PDF = "docs/ebf3/Ebf3_report.md", "docs/ebf3/Ebf3_report.pdf"
F = "figures/ebf3"


def fig(name, caption, width="100%"):
    p = f"{F}/{name}.png"
    if not os.path.exists(p):
        return f"*Figure pending: `{p}` not built yet.*\n"
    return f"![{caption}](../../{p}){{width={width}}}\n"


def table(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(("–" if pd.isna(v) else f"{v:.2f}") if isinstance(v, float)
                                     else ("–" if str(v) in ("", "nan") else str(v)) for v in r) + " |")
    return "\n".join(out) + "\n"


halves = json.load(open("results/ebf3/ep_oe_experimental_halves.json"))
ev = pd.read_csv("results/ebf3/regulatory_evidence.tsv", sep="\t")
ev_in = ev[ev["in_sim_window"].astype(str) == "True"]
ctcf = json.load(open("results/ebf3/ctcf_barriers.json"))
comp_path = "results/ebf3/model_comparison.csv"
comp = pd.read_csv(comp_path) if os.path.exists(comp_path) else None
sims_done = comp is not None and (comp["model"] != "experiment").any()

md = [f"""---
title: "Ebf3 / hs737: is loop extrusion enough?"
subtitle: "Experimental Hi-C and polymer simulations of chr7:135.8–137.7 Mb (mm10)"
date: "{datetime.date.today().isoformat()}"
geometry: margin=1.6cm
fontsize: 10pt
header-includes:
  - \\usepackage{{float}}
  - \\floatplacement{{figure}}{{H}}
mainfont: DejaVu Sans
colorlinks: true
---

# Question

VISTA enhancer **hs737** (chr7:136,083,204–136,084,349; mouse ortholog Rr169617,
whose knockout lowers Ebf3 ~20 %) lies **1.23 Mb** from the **Ebf3 promoter**
(TSS chr7:137,314,445, minus strand). Does cohesin loop extrusion with CTCF
barriers reproduce their contact frequency, or does an attractive
enhancer–promoter (E–P) interaction have to be added?

Statements are labelled as **observation** (measured), **prediction** (model
output), or **hypothesis**. A better fit from adding E–P attraction shows the
attraction is *sufficient* in the model, not that it is *necessary* in cells.
Hi-C is a population average; every result is an ensemble contact frequency.

# Data and their limits

| Data | Source | Cells | Role |
|---|---|---|---|
| Hi-C | GEO GSE161259 (HiCUP reanalysis of Bonev 2017) | ES; **neocortex-sorted** NPC, CN (`ncx_`) | experiment |
| CTCF, H3K27ac, H3K4me1/3, Input | GEO GSE96107 (Bonev 2017) | ES; **in vitro** NPC, CN | barriers, enhancer evidence |
| ATAC-seq | ENCODE ENCSR810HQR | bulk E14.5 forebrain | **proxy** for accessibility |
| Enhancers | VISTA (UCSC mm10) | — | annotation |

**Cell-source mismatch.** The Hi-C for NPC / CN is from neocortex, but the
ChIP is from in vitro cultures. Ebf3 is 30–79 FPKM in in vitro CN and only
1–4 FPKM in neocortex CN. Matched in vitro Hi-C must be reprocessed from SRA
(87 runs, ~786 GB; scoped, awaiting approval). Until then, every comparison of
the models (built on in vitro CN CTCF) with ncx_CN Hi-C carries this mismatch.

# Region

Simulated and displayed window: **chr7:135,800,000–137,700,000** (1.9 Mb,
1,900 × 1 kb monomers): 283 kb left of hs737, 386 kb right of the Ebf3 TSS,
enclosing the conserved CTCF clusters on both sides. The E–P contact is
measured at **10 kb**, the finest resolution with ≥ 30 reads in the E–P pixel
in every state (rule fixed before any data were looked at). The left edge sits
20–80 kb from a strong insulation boundary just outside the window
(docs/ebf3/01_hic_qc.md).

# Experimental Hi-C (observation)

{fig("fig1_experimental_hic", "Experimental Hi-C, all available Bonev states")}

{fig("fig1b_experimental_hic_zoom", "Ebf3 side at 5 kb", "90%")}

E–P observed/expected at 10 kb (expected = mean of the same diagonal inside the
window), with the two random halves of each dataset:

| State | O/E | half A | half B |
|---|---|---|---|
""" + "".join(f"| {s} | {h['full']:.2f} | {h['half_A']:.2f} | {h['half_B']:.2f} |\n" for s, h in halves.items()) + """
**Observation.** The hs737 × Ebf3 contact is enriched over distance-matched
background in both neural states (2.0–2.2×) and less in ES (1.35×), a larger
difference than between halves of one dataset. The strongest dot in the same
column is not at hs737 but at the CTCF pair 136,489(+) × 137,321(−).

# CTCF barriers

""", fig("03_ctcf_orientation", "Oriented CTCF sites"), f"""
FIMO (MA0139.1, p < 1e-4) oriented ES {ctcf['ES']['n']}, NPC {ctcf['NPC']['n']}
and CN {ctcf['CN']['n']} sites in the window. In CN, a − site 7 kb from the Ebf3
TSS and + sites flanking hs737 form convergent pairs that nearly bracket the E–P
pair, so extrusion alone may raise contacts near the E–P pixel. CN calls 4.3×
fewer CTCF peaks genome-wide than ES, most likely a peak-calling difference.

# Enhancer evidence

""", fig("04_regulatory_landscape", "Regulatory landscape"), """
""", table(ev_in[["element", "tier", "CN_H3K27ac_over_input", "CN_H3K4me1_over_input",
                   "ATAC_peak_foldenrich", "ATAC_max_foldchange"]].rename(columns={
    "CN_H3K27ac_over_input": "CN K27ac/in", "CN_H3K4me1_over_input": "CN K4me1/in",
    "ATAC_peak_foldenrich": "ATAC peak FE*", "ATAC_max_foldchange": "ATAC max FC*"})), """
\\* ATAC = **proxy** (bulk E14.5 forebrain, not these cells). hs737 has no
H3K27ac IDR peak in NPC or CN but has the strongest ATAC peak among VISTA
elements in the window (fold 9.9) and the highest CN H3K27ac/input ratio.

# Models

All share one polymer, one window, identical sampling (3 replicates × 10,000
frames × 28 tiled copies) and CN (in vitro) CTCF.

| Model | Cohesin + CTCF | E–P attraction | Force | Compare with |
|---|---|---|---|---|
| A | no | no | default | B |
| A0 | no | no | attractive, 0 kT | C, D |
| B | yes | no | default | A |
| C | no | 3 kT | attractive | A0 |
| D | yes | 3 kT | attractive | A0 |

Cohesin: lifetime 75,000, separation 240, CTCF capture 0.5 (polysim defaults,
not fitted). E–P attraction acts only between hs737 and the Ebf3 promoter.

# Simulations vs experiment (prediction)

"""]
if sims_done:
    e = comp[comp.metric == "ep_oe_10kb"].copy()
    e = e[["model", "state", "value", "ci_low", "ci_high", "n_replicates", "kind"]]
    md += [fig("fig2_model_panels", "Model maps vs ncx_CN Hi-C"),
           fig("fig3_ep_endpoint", "E–P endpoint per model", "70%"),
           "\nE–P O/E at 10 kb (models: mean, bootstrap 95 % CI over replicates; experiment: split-half range):\n\n",
           table(e)]
else:
    md += ["*Pending: the production simulations are running. This section fills in when "
           "`scripts/ebf3/plot_model_panels.py` has run.*\n"]

md += ["""
# What can and cannot be concluded

- If B reaches the experimental E–P O/E, extrusion with CTCF is **sufficient** in
  this model; it does not exclude an E–P attraction in cells.
- If only C or D reaches it, an E–P attraction is **sufficient** in the model to
  reproduce the measurement; it is not shown to be **necessary**, nor unique.
- Different parameter sets can give indistinguishable maps; the sensitivity
  sweeps (cohesin lifetime / density / capture, attraction energy, enhancer set)
  are needed before quoting a range.
- The neural Hi-C is neocortex; the model's CTCF and enhancer evidence are in
  vitro. Matched in vitro Hi-C is the decisive comparison.

# Open decisions

1. Tier-2 enhancer set: rule-computed (mm1590, hs737, hs699) vs listed (hs737,
   hs699, hs232). ATAC supports hs232 (peak, FE 5.5) over mm1590 (no peak).
2. Reprocess in vitro NPC / CN Hi-C from SRA (~786 GB, 1–2 days).
3. Left window edge: keep, or move to 135,700,000 to enclose the nearby boundary.
4. Sensitivity sweeps (cost: 12 GPU jobs per added condition).

*Sources, scripts and full methods: docs/ebf3/00–04, branch `ebf3-locus`.*
"""]

with open(OUT_MD, "w") as f:
    f.write("".join(md))
r = subprocess.run(["pandoc", OUT_MD, "-o", OUT_PDF, "--pdf-engine=xelatex",
                    "--resource-path=docs/ebf3", "-V", "tables=true"], capture_output=True, text=True)
print(r.stderr[-2000:] if r.returncode else f"wrote {OUT_PDF}")
