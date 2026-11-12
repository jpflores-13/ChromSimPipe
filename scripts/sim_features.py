"""
Turn per-condition options into the per-monomer inputs polysim takes.

Every option reads peaks (narrowPeak / BED) inside the active locus, maps them
to the locus' monomers, and repeats them in every copy of the tiled chain:

  ctcf_capture {min, max}  per-site stall probabilities scaled by CTCF peak signal
  sticky                   monomer types + an attraction matrix (polysim's
                           monomer_types / interaction_matrix)
  loading                  more cohesin loading at peaks (polysim's per-monomer
                           birth array)
  pausing                  slower cohesin stepping at peaks (polysim's per-monomer
                           pause array)

polysim's run() builds its cohesin arrays with polysim.extrusion.build_lef_arrays
and fills loading and pausing uniformly; `lef_profiles` wraps that function for
one run so the arrays carry the profiles. Nothing here changes polysim itself.
"""

import contextlib
from itertools import combinations_with_replacement

import numpy as np


def read_peaks(path, chrom, start, end, resolution, span="summit"):
    """Local monomer indices of the peaks on `chrom` inside [start, end).

    span "summit": the monomer holding each peak's summit (narrowPeak column 10,
    else the peak middle); "peak": every monomer the peak overlaps. Returns
    (sorted unique monomers, {monomer: max signal}); signal is narrowPeak column 7
    (signalValue), else column 5.
    """
    n = (end - start) // resolution
    signal = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            r = line.rstrip("\n").split("\t")
            if len(r) < 3 or r[0] != chrom:
                continue
            s, e = int(r[1]), int(r[2])
            if e <= start or s >= end:
                continue
            try:
                value = float(r[6]) if len(r) > 6 else float(r[4]) if len(r) > 4 else 1.0
            except ValueError:
                value = 1.0
            if span == "peak":
                first = max(0, (s - start) // resolution)
                last = min(n - 1, (e - 1 - start) // resolution)
                monomers = range(first, last + 1)
            else:
                try:
                    summit = s + int(r[9]) if len(r) > 9 and int(r[9]) >= 0 else (s + e) // 2
                except ValueError:
                    summit = (s + e) // 2
                m = (summit - start) // resolution
                monomers = [m] if 0 <= m < n else []
            for m in monomers:
                signal[m] = max(signal.get(m, value), value)
    return np.array(sorted(signal), dtype=int), signal


def tile(local, tiling):
    """Repeat local monomer indices in every copy of the tiled chain."""
    local = np.asarray(local, dtype=int)
    offsets = np.arange(tiling["n_tiles"]) * tiling["tile_size"] + tiling["padding"]
    return (offsets[:, None] + local[None, :]).ravel()


def capture_by_signal(signals, lo, hi):
    """Stall probability per site from its peak signal, by rank: the weakest site
    gets `lo`, the strongest `hi`, the rest evenly in between (ties share a rank).
    Ranks, not raw values, so the scale of the peak caller's signal doesn't matter."""
    signals = np.asarray(signals, dtype=float)
    if signals.size == 0:
        return signals
    if signals.size == 1 or np.all(signals == signals[0]):
        return np.full(signals.size, (lo + hi) / 2.0)
    order = signals.argsort(kind="stable")
    ranks = np.empty(signals.size)
    ranks[order] = np.arange(signals.size)
    for v in np.unique(signals):            # average ranks of ties
        tie = signals == v
        ranks[tie] = ranks[tie].mean()
    return lo + (hi - lo) * ranks / (signals.size - 1)


def stall_dicts(positions, orientations, probs, groups=None, combine="max"):
    """polysim's ctcf_left / ctcf_right {monomer: probability} for sites of one orientation
    sharing a monomer:
      combine="max"          keep the highest probability
      combine="independent"  1 - prod(1 - p) over distinct `groups` (peak names): the chance
                             at least one bound CTCF stalls the leg (Harvey's supplement,
                             8.3); sites of the same group (motifs of one peak) count once,
                             at their highest probability."""
    if combine not in ("max", "independent"):
        raise ValueError(f"combine must be 'max' or 'independent', not {combine!r}")
    positions = np.asarray(positions, int)
    if groups is None:
        groups = range(len(positions))
    best = {}                                   # (side, monomer, group) -> probability
    for m, o, p, g in zip(positions, orientations, probs, groups):
        k = (+1 if o == +1 else -1, int(m), g)
        best[k] = max(best.get(k, 0.0), float(p))
    left, right = {}, {}
    for (o, m, _), p in best.items():
        side = left if o == +1 else right
        if combine == "max":
            side[m] = max(side.get(m, 0.0), p)
        else:
            side[m] = 1.0 - (1.0 - side.get(m, 0.0)) * (1.0 - p)
    return left, right


def sticky_types(npoly, positions, energies, cross=()):
    """Monomer types and attraction matrix for sticky factors.

    positions: {factor: global monomer indices}; energies: {factor: kT between two
    of its own sites}; cross: [{factors: [A, B], energy}]. A monomer bound by
    several factors gets a type of its own for that combination; two types attract
    with the strongest energy of any factor pair they hold. Type 0 is unmarked
    chromatin. Returns (types array, matrix, labels).
    """
    combo = [frozenset() for _ in range(npoly)]
    for factor, pos in positions.items():
        for m in np.asarray(pos, int):
            combo[m] = combo[m] | {factor}
    labels = [frozenset()] + sorted({c for c in combo if c}, key=lambda c: sorted(c))
    index = {c: i for i, c in enumerate(labels)}
    types = np.array([index[c] for c in combo], dtype=int)

    pair = {(f, f): float(e) for f, e in energies.items()}
    for x in cross:
        a, b = x["factors"]
        pair[(a, b)] = pair[(b, a)] = float(x["energy"])
    k = len(labels)
    matrix = np.zeros((k, k))
    for i, j in combinations_with_replacement(range(k), 2):
        e = max((pair.get((a, b), 0.0) for a in labels[i] for b in labels[j]), default=0.0)
        matrix[i, j] = matrix[j, i] = e
    return types, matrix, ["+".join(sorted(c)) or "none" for c in labels]


@contextlib.contextmanager
def lef_profiles(loading=None, pausing=None):
    """Within the block, polysim builds its cohesin arrays with local loading and
    pausing: loading = (global monomers, boost) multiplies the loading weight there;
    pausing = (global monomers, slowdown) divides the stepping probability there."""
    from polysim import extrusion
    original = extrusion.build_lef_arrays
    if not loading and not pausing:
        yield
        return

    def patched(length, *args, **kwargs):
        arrays = original(length, *args, **kwargs)
        if loading:
            pos, boost = loading
            arrays["birth"][np.asarray(pos, int)] *= boost
        if pausing:
            pos, slowdown = pausing
            step = 1.0 - arrays["pause"]
            arrays["pause"][np.asarray(pos, int)] = 1.0 - step[np.asarray(pos, int)] / slowdown
        return arrays

    extrusion.build_lef_arrays = patched
    try:
        yield
    finally:
        extrusion.build_lef_arrays = original
