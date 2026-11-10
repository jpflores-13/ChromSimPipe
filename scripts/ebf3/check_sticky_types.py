"""Stage 4 check: the monomer types / interaction matrix / labels that
run_simulation_shard.py passes to polysim for each sticky condition of a config,
built exactly as that script does (sim_features.read_peaks / tile / sticky_types).
Usage: CHROMSIM_CONFIG=config/ChromSimConfig.ebf3.yaml python scripts/ebf3/check_sticky_types.py"""
import os, sys
sys.path.insert(0, os.getcwd())
import numpy as np
import configs.parameters as P
from scripts import sim_features as feat

for cond in P.SIMULATION_CONDITIONS:
    prm = cond["params"]
    if not prm.get("sticky"):
        continue
    files = P.condition_peak_files(cond)
    positions, energies, local = {}, {}, {}
    for f in prm["sticky"]:
        loc, _ = feat.read_peaks(files[f"sticky:{f['factor']}"], P.CHROM, P.REGION_START,
                                 P.REGION_END, P.RESOLUTION, f["span"])
        local[f["factor"]] = loc.tolist()
        positions[f["factor"]] = feat.tile(loc, P.TILING)
        energies[f["factor"]] = f["energy"]
    types, matrix, labels = feat.sticky_types(int(P.TILING["chrom_size"]), positions, energies,
                                              prm.get("sticky_cross", []))
    print(f"== {cond['name']}  (chain {P.TILING['chrom_size']} monomers, {P.TILING['n_tiles']} copies)")
    for k, v in local.items():
        print(f"  {k}: local monomers {v} -> bp {[P.REGION_START + m * P.RESOLUTION for m in v]}")
    shared = set.intersection(*(set(v) for v in local.values())) if len(local) > 1 else set()
    print(f"  E and P summits sharing a bin: {sorted(shared) or 'none'}")
    print(f"  labels: {labels}")
    print(f"  matrix (kT):\n{np.array2string(matrix, precision=2)}")
    u, n = np.unique(types, return_counts=True)
    print(f"  type counts on the chain: {dict(zip(u.tolist(), n.tolist()))}")
    print(f"  types in copy 0 at those monomers: "
          f"{ {m: int(types[P.TILING['padding'] + m]) for v in local.values() for m in v} }")
