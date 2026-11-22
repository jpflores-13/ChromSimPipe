#!/usr/bin/env python
"""
Export the Hi-C and simulated maps behind the cell-type x model rectangle figure
(scripts/ebf3/plot_celltype_rectangles.R) as sparse upper-triangle tables.

  python scripts/ebf3/export_panel_matrices.py <locus> [variant]

<locus> names results/<locus>/panel_locus.json (ebf3, bdnf, bdnf_wide, ...). [variant]
(e.g. occ, tier) takes the B and D columns from B_extrusion_<variant> / D_combined_t1_<variant>
(per-site CTCF stall probabilities) and writes to results/<locus>/panels_<variant>/.

Locus definition: results/<locus>/panel_locus.json
  chrom, sim [start, end] (the simulated region), analysis (pipeline analysis folder),
  hic {ES, NPC, CN: mcool}, hic_labels {cell: label}, gene [name, start, end, strand],
  enhancer [name, start, end] | null, ep [enhancer_bp, promoter_bp] | null,
  anchors [[start, end], ...] (highlight bands), signals {track: {cell: bigWig | null}}
  (drawn by the R script), view [start, end] (optional).

View: unless `view` is given, the widest window inside the simulated region whose
plotgardener rectangle (height/width = RATIO) is filled by simulated data at its top
corners (they need RATIO x span of data beyond each edge), centred on gene + enhancer.

For each cell type (ES, NPC, CN) x column (Hi-C, A, A0, B, C, D):
  results/<locus>/panels/<cell>_<column>.tsv   pos1  pos2  value (10 kb)
and, when `ep` is set, results/<locus>/panels/loops.bedpe (the E-P pixel, for annoPixels);
when `arches` {promoter: [s, e], peaks: narrowPeak, merge_bp, label} is set,
results/<locus>/panels/arches.bedpe (promoter x merged enhancer clusters, for plotPairsArches).
value = contact / the map's own mean contact at 50 kb (sim and Hi-C units differ), so
every panel shares one colour scale. Sim maps are the pooled 3-replicate maps from
analyzeCondition; a map that does not exist yet is listed as pending.
E-P O/E: the 10 kb E x P pixel / mean of its diagonal over oe_window (Ebf3: the pre-registered
chr7:135.8-137.7 Mb, scripts/ebf3/hic_qc.py), identical for Hi-C and every model.
Run from the repo root with the ebf3_analysis env.
"""
import json, os, sys
import numpy as np
import cooler

BIN, RES = 10_000, 1000
F = BIN // RES
RATIO = 1.0 / 2.3                     # panel height / width in the R script (H / W)
MODELS = [("A_polymer", "polymer"), ("A0_polymer_attr", "polymer_attr"), ("B_extrusion", "extrusion"),
          ("C_sticky_t1", "sticky_t1"), ("D_combined_t1", "combined_t1")]
CELLS = ("ES", "NPC", "CN")


def coarsen(m):
    n = m.shape[0] // F
    m = np.nan_to_num(m[: n * F, : n * F])
    return m.reshape(n, F, n, F).sum(axis=(1, 3))


def choose_view(L):
    if L.get("view"):
        return tuple(L["view"])
    s0, s1 = L["sim"]
    feats = [L["gene"][1], L["gene"][2]] + ([L["enhancer"][1], L["enhancer"][2]] if L.get("enhancer") else [])
    lo, hi = min(feats), max(feats)
    span = (s1 - s0) / (1 + 2 * RATIO) // BIN * BIN
    m = RATIO * span
    a = min(max((lo + hi) / 2 - span / 2, s0 + m), s1 - m - span)
    a = -(-a // BIN) * BIN
    b = a + span
    if lo < a or hi > b:
        sys.exit(f"gene/enhancer ({lo}-{hi}) do not fit a filled view inside the simulated region")
    if L.get("ep") and abs(L["ep"][1] - L["ep"][0]) > 2 * RATIO * span:
        sys.exit("the E-P pixel would sit above the rectangle")
    return int(a), int(b)


def main():
    name = sys.argv[1]
    variant = sys.argv[2] if len(sys.argv) > 2 else ""
    L = json.load(open(f"results/{name}/panel_locus.json"))
    models = [(c, p) if c not in ("B_extrusion", "D_combined_t1") or not variant
              else (f"{c}_{variant}", f"{p}_{variant}") for c, p in MODELS]
    chrom, (s0, s1) = L["chrom"], L["sim"]
    ws, we = choose_view(L)
    out = f"results/{name}/panels" + (f"_{variant}" if variant else "")
    os.makedirs(out, exist_ok=True)
    margin = int(RATIO * (we - ws)) + BIN
    h0, h1 = max(0, ws - margin) // BIN * BIN, -(-(we + margin) // BIN) * BIN   # Hi-C export range

    panels, vals, depth = [], [], {}
    for cell in CELLS:
        c = cooler.Cooler(f"{L['hic'][cell]}::resolutions/{BIN}")
        depth[cell] = int(c.info.get("sum") or c.pixels()[:]["count"].sum())
        hic_ref = c.matrix(balance=True).fetch(f"{chrom}:{s0}-{s1}")     # 50 kb / O/E references
        hic_out = c.matrix(balance=True).fetch(f"{chrom}:{h0}-{h1}")
        cols = [("HiC", hic_ref, hic_out, h0)]
        for (cond, params), (col, _) in zip(models, MODELS):
            cname = cond if cell == "CN" else f"{cond}_{cell}"
            p = f"{L['analysis']}/{cname}/{params}_ctcf-{cell}_30000blk_pooled_contact_map.npy"
            m = coarsen(np.load(p)) if os.path.exists(p) else None
            cols.append((col, m, m, s0))
        for col, ref, m, origin in cols:
            rec = dict(cell=cell, column=col, file=None, ep_oe=None)
            if m is not None:
                ref = np.where(np.isfinite(ref), ref, np.nan)
                m = np.where(np.isfinite(m), m, np.nan)
                d50 = np.nanmean(np.diagonal(ref, 50_000 // BIN))
                if L.get("ep"):
                    o0, o1 = L.get("oe_window", L["sim"])     # diagonal-mean window of the endpoint
                    r = ref[(o0 - s0) // BIN:(o1 - s0) // BIN, (o0 - s0) // BIN:(o1 - s0) // BIN]
                    i, j = (L["ep"][0] - o0) // BIN, (L["ep"][1] - o0) // BIN
                    rec["ep_oe"] = float(r[i, j] / np.nanmean(np.diagonal(r, j - i)))
                v = m / d50
                iu, ju = np.triu_indices(v.shape[0])
                keep = np.isfinite(v[iu, ju]) & (v[iu, ju] > 0)
                f = f"{out}/{cell}_{col}.tsv"
                np.savetxt(f, np.c_[origin + iu[keep] * BIN, origin + ju[keep] * BIN, v[iu, ju][keep]],
                           fmt=["%d", "%d", "%.6g"], delimiter="\t", header="pos1\tpos2\tvalue", comments="")
                rec["file"] = f
                a, b = (ws - origin) // BIN, (we - origin) // BIN
                w = v[a:b, a:b][np.triu_indices(b - a, k=1)]     # off-diagonal sets the scale
                vals.append(w[np.isfinite(w) & (w > 0)])
            panels.append(rec)

    loops = None
    if L.get("ep"):   # the E-P pixel as a bedpe (10 kb bins), read by annoPixels
        e, pr = sorted(L["ep"])
        loops = f"{out}/loops.bedpe"
        en = L["enhancer"][0] if L.get("enhancer") else "E"
        with open(loops, "w") as fh:
            fh.write(f"{chrom}\t{e // BIN * BIN}\t{e // BIN * BIN + BIN}\t"
                     f"{chrom}\t{pr // BIN * BIN}\t{pr // BIN * BIN + BIN}\t{en}-{L['gene'][0]}_promoter\n")

    arches = None
    if L.get("arches"):   # promoter x candidate-enhancer clusters, for plotPairsArches
        A = L["arches"]
        p0, p1 = A["promoter"]
        arches = f"{out}/arches.bedpe"

        def cluster(peaks_file):          # merge peaks closer than merge_bp
            peaks = sorted((int(r[1]), int(r[2]), r[3], float(r[6])) for r in
                           (l.split("\t") for l in open(peaks_file)) if r[0] == chrom)
            cl = []
            for st, en, nm, sig in peaks:
                if cl and st - cl[-1]["end"] <= A.get("merge_bp", 20_000):
                    c = cl[-1]; c["end"] = max(c["end"], en); c["names"].append(nm)
                    c["h3k27ac"] |= "H3K27ac" in nm; c["max"] = max(c["max"], sig)
                else:
                    cl.append(dict(start=st, end=en, names=[nm], h3k27ac="H3K27ac" in nm, max=sig))
            return cl

        ## `peaks` is one narrowPeak (legacy: evidence from the peak names) or {evidence: narrowPeak},
        ## one evidence type per file, clustered separately (no union)
        sets = ([(None, cluster(A["peaks"]))] if isinstance(A["peaks"], str)
                else [(ev, cluster(f)) for ev, f in A["peaks"].items()])
        with open(arches, "w") as fh:
            fh.write("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tevidence\tn_peaks\tmax_signal\n")
            for ev0, clusters in sets:
                for c in clusters:
                    a, b = sorted([(p0, p1), (c["start"], c["end"])])
                    ev = ev0 or ("CN_H3K27ac" if c["h3k27ac"] else "ATAC_only")
                    nm = (f"{c['names'][0].split('_')[1]}-{c['names'][-1].split('_')[1]}" if ev0 is None
                          else f"{c['names'][0].split('_')[-1]}-{c['names'][-1].split('_')[-1]}")
                    fh.write(f"{chrom}\t{a[0]}\t{a[1]}\t{chrom}\t{b[0]}\t{b[1]}\t"
                             f"{nm}\t{ev}\t{len(c['names'])}\t{c['max']:.1f}\n")

    allv = np.concatenate(vals)
    meta = dict(locus=name, variant=variant, chrom=chrom, start=ws, end=we, sim=[s0, s1], bin=BIN, ratio=RATIO,
                assembly=L.get("assembly", "mm10"), depth=depth, hic_labels=L.get("hic_labels", {}),
                gene=dict(zip(("name", "start", "end", "strand"), L["gene"])),
                enhancer=dict(zip(("name", "start", "end"), L["enhancer"])) if L.get("enhancer") else None,
                ep=L.get("ep"), loops=loops, arches=arches,
                arches_label=(L.get("arches") or {}).get("label"), anchors=L.get("anchors", []), signals=L.get("signals", {}),
                signal_labels=L.get("signal_labels", {}),
                zrange=[float(np.percentile(allv, 10)), float(np.percentile(allv, 97))],
                panels=panels)
    json.dump(meta, open(f"{out}/meta.json", "w"), indent=1)
    done = sum(p["file"] is not None for p in panels)
    print(f"{name}: view {chrom}:{ws:,}-{we:,} ({(we - ws) / 1e6:.2f} Mb) in sim {s0:,}-{s1:,}; "
          f"{done}/{len(panels)} panels -> {out}")


if __name__ == "__main__":
    main()
