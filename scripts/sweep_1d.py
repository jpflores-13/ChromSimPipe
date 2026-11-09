#!/usr/bin/env python
"""
Cheap 1D parameter sweep: cohesin on the CTCF track only, no 3D polymer, no GPU.

For one locus and one condition from the config, runs polysim's 1D cohesin
engine for every combination of the values given (lifetime, separation,
ctcf_capture, vlef, life_boost_stalled) and reports, for each:

  looped fraction   per loop: the fraction of time a cohesin holds its two
                    anchors (left leg at one anchor, right leg at the other,
                    within --tol monomers), pooled over all copies of the locus
  mean loop size    average distance between the legs of a bound cohesin, kb
  coverage          fraction of the chain inside at least one cohesin loop

It is a screen: which settings could give loops at all, how long they last,
and how a perturbation shifts them, before spending GPU time on 3D runs. The
looped fraction is the 1D quantity live-imaging studies measure for one locus
(e.g. Gabriele et al. 2022). Sticky sites act in 3D only and are ignored here;
loading and pausing profiles are applied.

Run with the simulation environment's Python (see sweep_1d.sh):
    bash sweep_1d.sh --locus chr1_fig1 --condition cond1_baseline \\
        --lifetime 37500,75000,150000 --ctcf-capture 0.25,0.5 --out sweep.tsv
Loops default to every convergent CTCF pair (+ ... -) in the locus between
--min-kb and --max-kb; --loops takes a TSV with chrom, anchor1, anchor2 (bp).
"""

import argparse
import itertools
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)


def floats(text):
    return [float(x) for x in str(text).split(",") if x.strip()]


def score_frame(L, R, a_all, b_all, tol, npoly):
    """One snapshot of cohesin legs L, R (monomers). a_all, b_all: (loops, copies)
    anchor monomers. Returns (fraction of copies looped per loop, mean loop size,
    fraction of the chain inside a loop)."""
    L, R = np.asarray(L, int), np.asarray(R, int)
    near_a = np.abs(L[:, None, None] - a_all[None]) <= tol
    near_b = np.abs(R[:, None, None] - b_all[None]) <= tol
    looped = (near_a & near_b).any(axis=0).mean(axis=1)
    cover = np.zeros(npoly + 1, dtype=int)
    np.add.at(cover, L, 1)
    np.add.at(cover, R, -1)
    return looped, float(np.mean(R - L)), float(np.mean(np.cumsum(cover)[:npoly] > 0))


def label(v):
    return "{min:g}-{max:g} by signal".format(**v) if isinstance(v, dict) else f"{v:g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", default=None, help="Locus name from the config (default: the first)")
    ap.add_argument("--condition", required=True, help="Condition giving the CTCF sites and defaults")
    ap.add_argument("--lifetime", default=None, help="Comma-separated values (default: the condition's)")
    ap.add_argument("--separation", default=None)
    ap.add_argument("--ctcf-capture", default=None)
    ap.add_argument("--vlef", default=None)
    ap.add_argument("--life-boost-stalled", default=None)
    ap.add_argument("--loops", default=None, help="TSV with chrom, anchor1, anchor2 columns (bp)")
    ap.add_argument("--min-kb", type=float, default=50)
    ap.add_argument("--max-kb", type=float, default=1000)
    ap.add_argument("--tol", type=int, default=2, help="Monomers a leg may sit from an anchor")
    ap.add_argument("--frames", type=int, default=2000, help="Snapshots recorded per setting")
    ap.add_argument("--every", type=int, default=None,
                    help="Cohesin steps between snapshots (default: lifetime / 20)")
    ap.add_argument("--lifetimes-equil", type=float, default=10,
                    help="Equilibration length, in (boosted) cohesin lifetimes")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    if args.locus:
        os.environ["CHROMSIM_LOCUS"] = args.locus
    import configs.parameters as p
    import sim_features as feat
    from polysim import extrusion
    import run_simulation_shard as rs

    cond = p.get_condition(args.condition)
    base = cond["params"]
    if base["separation"] is None:
        sys.exit(f"{args.condition} has no cohesin (separation: none); nothing to sweep in 1D")
    loc = p.LOCI[p.ACTIVE_LOCUS]
    pos, ori = p.get_ctcf_arrays(cond["ctcf_type"])
    tpos, tori = p.get_tiled_ctcf_arrays(cond["ctcf_type"])
    npoly, tiles = int(p.TILING["chrom_size"]), p.TILING
    offsets = np.arange(tiles["n_tiles"]) * tiles["tile_size"] + tiles["padding"]

    if args.loops:
        loops = []
        with open(args.loops) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(header, line.rstrip("\n").split("\t")))
                if r.get("chrom") != loc["chrom"]:
                    continue
                a = (int(r["anchor1"]) - loc["region_start"]) // p.RESOLUTION
                b = (int(r["anchor2"]) - loc["region_start"]) // p.RESOLUTION
                if 0 <= a < b < p.N_MONOMERS:
                    loops.append((a, b))
    else:
        lo, hi = args.min_kb * 1000 / p.RESOLUTION, args.max_kb * 1000 / p.RESOLUTION
        loops = [(int(a), int(b)) for (a, oa), (b, ob) in itertools.combinations(zip(pos, ori), 2)
                 if oa == +1 and ob == -1 and lo <= b - a <= hi]
    if not loops:
        sys.exit("No loops to score: no convergent CTCF pairs in range, or none from --loops")

    features, _ = rs.condition_features(cond)
    grid = dict(
        lifetime=floats(args.lifetime) if args.lifetime else [base["lifetime"]],
        separation=floats(args.separation) if args.separation else [base["separation"]],
        ctcf_capture=floats(args.ctcf_capture) if args.ctcf_capture else [base["ctcf_capture"]],
        vlef=floats(args.vlef) if args.vlef else [base.get("vlef", p.POLYSIM["vlef"])],
        life_boost_stalled=floats(args.life_boost_stalled) if args.life_boost_stalled
        else [base.get("life_boost_stalled", p.POLYSIM["life_boost_stalled"])],
    )
    print(f"{p.ACTIVE_LOCUS}: {len(loops)} loops, {tiles['n_tiles']} copies, "
          f"{np.prod([len(v) for v in grid.values()])} settings")

    rows = []
    for values in itertools.product(*grid.values()):
        s = dict(zip(grid, values))
        params = {**base, "ctcf_capture": s["ctcf_capture"]}
        left, right = rs.ctcf_stall_sites(params, cond["ctcf_type"], tpos, tori)
        stall_left, stall_right = extrusion.build_stall_arrays(npoly, ctcf_left=left,
                                                               ctcf_right=right, stall_prob=None)
        rs.seed_all(args.seed)
        with feat.lef_profiles(features.get("loading"), features.get("pausing")):
            arrays = extrusion.build_lef_arrays(npoly, lifetime=s["lifetime"], vlef=s["vlef"],
                                                stall_left=stall_left, stall_right=stall_right,
                                                life_boost_stalled=s["life_boost_stalled"])
        n_lefs = int(npoly // s["separation"])
        smc = extrusion.make_translocator(arrays, n_lefs)
        t0 = time.time()
        smc.steps(int(args.lifetimes_equil * s["lifetime"] * s["life_boost_stalled"]))
        every = args.every or max(1, int(s["lifetime"] / 20))
        looped = np.zeros(len(loops))
        sizes, covered = [], []
        a_all = np.array([a for a, _ in loops])[:, None] + offsets[None, :]
        b_all = np.array([b for _, b in loops])[:, None] + offsets[None, :]
        for _ in range(args.frames):
            smc.steps(every)
            L, R = smc.getLEFs()
            f, size, cov = score_frame(L, R, a_all, b_all, args.tol, npoly)
            looped += f
            sizes.append(size)
            covered.append(cov)
        looped /= args.frames
        row = {**{k: label(v) for k, v in s.items()}, "mean_loop_kb": np.mean(sizes) * p.RESOLUTION / 1000,
               "coverage": np.mean(covered), "seconds": time.time() - t0}
        for (a, b), f in zip(loops, looped):
            row[f"looped_{loc['region_start'] + a * p.RESOLUTION}_{loc['region_start'] + b * p.RESOLUTION}"] = f
        rows.append(row)
        print("  " + "  ".join(f"{k}={label(v)}" for k, v in s.items())
              + f"  mean loop {row['mean_loop_kb']:.0f} kb, coverage {row['coverage']:.2f}, "
              f"max looped fraction {looped.max():.3f}")

    cols = list(rows[0])
    with open(args.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
