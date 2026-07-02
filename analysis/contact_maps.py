#!/usr/bin/env python
"""
Generate contact maps from simulation trajectories and compare with experimental Hi-C.

Reads 3D polymer conformations (from polychrom or OpenMM) and computes:
  1. Simulated contact maps
  2. P(s) contact probability decay curves
  3. Insulation score profiles
  4. Comparison metrics with experimental Hi-C

Usage:
    python contact_maps.py --sim-dir ../results/mESC_rep0 --exp-hic ../data/hic_mESC_Sox2.npy
"""

import os
import sys
import argparse
import logging
import json

if os.environ.get("CONDA_DEFAULT_ENV", "") != "cohesin_sim":
    sys.stderr.write(
        f"[env] active conda env is '{os.environ.get('CONDA_DEFAULT_ENV') or 'none'}', "
        "expected 'cohesin_sim'. run: conda activate cohesin_sim\n"
    )
    sys.exit(1)

import numpy as np
# ---------------------------------------------------------------------------
# scipy compatibility shim.
#
# `cKDTree` and `KDTree` have been the same class since scipy 1.6 (2021),
# and scipy >= 1.14 finally dropped the `cKDTree` alias. Fall back to
# `KDTree` when the legacy name isn't available so the pipeline runs on
# any modern scipy without requiring an env reinstall.
# ---------------------------------------------------------------------------
try:
    from scipy.spatial import cKDTree  # scipy < 1.14 (and broken installs that happen to still expose it)
except ImportError:  # pragma: no cover — exercised only on scipy >= 1.14
    from scipy.spatial import KDTree as cKDTree
from typing import Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import N_MONOMERS, RESOLUTION, SIM_RUN, TILING

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =============================================================================
# CONTACT MAP COMPUTATION
# =============================================================================

def _find_contacts(args):
    """Process one frame: return sparse contact pairs. Module-level for pickling."""
    coords, contact_radius = args
    tree = cKDTree(coords)
    pairs = tree.query_pairs(contact_radius, output_type="ndarray")
    return pairs


def _block_file_tile_contact_worker(args):
    """Hansen-style worker: open one HDF5 block file, walk every frame × tile
    sub-region, return a partial (N_locus, N_locus) contact map and the
    sub-region count. Nothing crosses the parent heap except that small result.

    Mirrors the pattern in polychrom.contactmaps.averageContacts where each
    worker calls ``loadFunction(self.filenames[self.i])`` one block at a time
    rather than receiving pre-loaded coordinates.
    """
    import h5py
    import numpy as _np
    try:
        from scipy.spatial import cKDTree as _cKDTree
    except ImportError:
        from scipy.spatial import KDTree as _cKDTree

    block_file, tile_offsets, N_locus, contact_radius = args
    tile_offsets = _np.asarray(tile_offsets, dtype=int)
    local_map = _np.zeros((N_locus, N_locus), dtype=_np.float64)
    n_sub = 0
    with h5py.File(block_file, "r") as hf:
        keys = sorted(hf.keys(), key=lambda x: int(x) if x.isdigit() else 0)
        for key in keys:
            grp = hf[key]
            if not (isinstance(grp, h5py.Group) and "pos" in grp):
                continue
            frame = _np.asarray(grp["pos"], dtype=_np.float64)
            for off in tile_offsets:
                sub = frame[off : off + N_locus]
                if sub.shape[0] != N_locus:
                    continue
                tree = _cKDTree(sub)
                pairs = tree.query_pairs(contact_radius, output_type="ndarray")
                if len(pairs) > 0:
                    local_map[pairs[:, 0], pairs[:, 1]] += 1.0
                    local_map[pairs[:, 1], pairs[:, 0]] += 1.0
                n_sub += 1
            # Free per-frame handles before next key to keep worker heap flat.
            del frame
    return local_map, n_sub


def _extract_tiles_file_parallel(
    block_files,
    contact_radius: float,
    n_jobs: int,
    *,
    tile_offsets,
    N_locus: int,
) -> np.ndarray:
    """File-parallel contact-map accumulator.

    Each worker gets one block_file path and returns a (N_locus, N_locus)
    partial map. The parent only sums those partials — no trajectory data
    ever enters the parent heap.

    This is the structural fix recommended in the Hansen-lab post-mortem
    (2026-04-22): parent heap stays at ~tens of MB regardless of trajectory
    length, eliminating the fork-CoW amplification that was causing
    2.8M-tile runs to grow the parent to 100+ GB.
    """
    from multiprocessing import get_context
    import gc, ctypes

    tile_offsets_list = [int(x) for x in tile_offsets]
    tasks = [
        (bf, tile_offsets_list, int(N_locus), float(contact_radius))
        for bf in block_files
    ]
    logger.info(
        f"  File-parallel contact extraction: {len(block_files)} blocks × "
        f"{len(tile_offsets_list)} tiles, {n_jobs} workers"
    )

    final_map = np.zeros((N_locus, N_locus), dtype=np.float64)
    total_n = 0

    # Fork context is fine here: parent heap is near-empty at this point
    # (no (T, N, 3) trajectory held), so inherited CoW pages are tiny.
    ctx = get_context("fork")
    done = 0
    try:
        libc = ctypes.CDLL("libc.so.6")
    except OSError:
        libc = None

    with ctx.Pool(n_jobs) as pool:
        for local_map, local_n in pool.imap_unordered(
            _block_file_tile_contact_worker, tasks, chunksize=1,
        ):
            final_map += local_map
            total_n += local_n
            # Explicit drop: prevents any per-iteration retention in
            # multiprocessing's IMapUnorderedIterator chunk cache, which was
            # the proximate cause of the 2026-04-22 streaming-path OOM
            # (job 9380631 leaked 1 pairs ndarray per frame; same pattern
            # could in principle retain one 32 MB local_map per block).
            del local_map
            done += 1
            if done % 25 == 0 or done == len(tasks):
                logger.info(
                    f"  [file-parallel] {done}/{len(tasks)} blocks, "
                    f"{total_n:,} sub-regions accumulated"
                )
                gc.collect()
                if libc is not None:
                    libc.malloc_trim(0)
    if total_n == 0:
        raise RuntimeError("File-parallel extraction processed 0 sub-regions.")
    return final_map / total_n


def compute_contact_map_from_conformations(
    conformations,
    contact_radius: float = 3.0,
    n_jobs: int = 1,
    n_frames: int = None,
) -> np.ndarray:
    """
    Compute average contact map from an ensemble of 3D conformations.

    Accepts either a list or a lazy generator of (N, 3) coordinate arrays.
    Using a generator avoids holding all frames in memory simultaneously.

    Parameters
    ----------
    conformations : iterable of np.ndarray
        Each element is an (N, 3) array of monomer positions.
        Can be a list or a generator.
    contact_radius : float
        Distance threshold for contact detection.
    n_jobs : int
        Number of parallel workers (default 1 = serial).
    n_frames : int, optional
        Total frame count — used only for progress logging when conformations
        is a generator (its length cannot be known upfront).

    Returns
    -------
    contact_map : np.ndarray of shape (N, N)
        Average contact frequency matrix.
    """
    from multiprocessing import Pool

    contact_map = None
    frame_count = 0
    total_str = f"/{n_frames}" if n_frames is not None else ""

    def _task_iter(confs, radius):
        for coords in confs:
            yield (coords, radius)

    if n_jobs > 1:
        logger.info(f"  Processing frames with {n_jobs} workers...")
        with Pool(n_jobs) as pool:
            for pairs in pool.imap(_find_contacts, _task_iter(conformations, contact_radius),
                                   chunksize=100):
                if contact_map is None:
                    # First frame: infer N from pairs or wait for next frame
                    # We need N — peek ahead by tracking coords shape externally.
                    # N is injected via the wrapper below.
                    pass
                if pairs is not None and contact_map is not None and len(pairs) > 0:
                    contact_map[pairs[:, 0], pairs[:, 1]] += 1
                    contact_map[pairs[:, 1], pairs[:, 0]] += 1
                frame_count += 1
                if frame_count % 500 == 0:
                    logger.info(f"  Frame {frame_count}{total_str}")
    else:
        logger.info(f"  Processing frames (serial)...")
        for coords in conformations:
            if contact_map is None:
                contact_map = np.zeros((coords.shape[0], coords.shape[0]), dtype=np.float64)
            pairs = _find_contacts((coords, contact_radius))
            if len(pairs) > 0:
                contact_map[pairs[:, 0], pairs[:, 1]] += 1
                contact_map[pairs[:, 1], pairs[:, 0]] += 1
            frame_count += 1
            if frame_count % 500 == 0:
                logger.info(f"  Frame {frame_count}{total_str}")

    if contact_map is None or frame_count == 0:
        raise RuntimeError("No frames processed — conformations iterable was empty.")

    contact_map /= frame_count
    return contact_map


def _accumulate_contacts(N: int, contact_radius: float, frame_iter, n_jobs: int,
                         n_frames: int = None) -> np.ndarray:
    """
    Internal helper: accumulate a contact map from a frame iterator.

    Works with both lists and generators. Returns the unnormalized sum
    and the frame count so the caller can normalize.

    Pool hardening (to keep memory bounded on long trajectories):
      - ``maxtasksperchild=500`` recycles each worker after 500 tiles so
        glibc-malloc fragmentation / accumulated per-process state does not
        grow unbounded for multi-million-task runs (2.8M sub-regions for a
        100k-frame 28-tile trajectory).
      - ``imap_unordered`` (instead of imap) lets results arrive as soon as
        any worker finishes, so the main thread drains the result queue
        continuously. Order does not matter because we only accumulate
        symmetric +=1 updates into a contact map.
      - ``chunksize=20`` keeps at most ~320 tiles (16 workers × 20) in
        flight through the pickle/IPC queue at any moment.
    """
    from multiprocessing import Pool

    contact_map = np.zeros((N, N), dtype=np.float64)
    frame_count = 0
    total_str = f"/{n_frames}" if n_frames is not None else ""

    def task_iter():
        for coords in frame_iter:
            yield (coords, contact_radius)

    if n_jobs > 1:
        import gc, resource
        with Pool(n_jobs, maxtasksperchild=500) as pool:
            for pairs in pool.imap_unordered(_find_contacts, task_iter(), chunksize=20):
                if len(pairs) > 0:
                    contact_map[pairs[:, 0], pairs[:, 1]] += 1
                    contact_map[pairs[:, 1], pairs[:, 0]] += 1
                frame_count += 1
                if frame_count % 500 == 0:
                    logger.info(f"  Frame {frame_count}{total_str}")
                # Periodically force Python GC + log parent RSS. Cheap
                # (runs every ~10k frames, ~30s real time) and gives us
                # the one diagnostic we were missing: *is parent RSS
                # growing linearly with frames* (real leak) or plateauing
                # (expected steady state).
                if frame_count % 10_000 == 0:
                    gc.collect()
                    parent_rss_mb = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss / 1024
                    child_rss_mb = resource.getrusage(
                        resource.RUSAGE_CHILDREN).ru_maxrss / 1024
                    logger.info(f"  [mem] frame {frame_count}{total_str}: "
                                f"parent peak RSS = {parent_rss_mb:.0f} MB, "
                                f"children peak RSS = {child_rss_mb:.0f} MB")
    else:
        for coords in frame_iter:
            pairs = _find_contacts((coords, contact_radius))
            if len(pairs) > 0:
                contact_map[pairs[:, 0], pairs[:, 1]] += 1
                contact_map[pairs[:, 1], pairs[:, 0]] += 1
            frame_count += 1
            if frame_count % 500 == 0:
                logger.info(f"  Frame {frame_count}{total_str}")

    return contact_map, frame_count


def extract_tiles_and_average(
    conformations,
    contact_radius: float = 3.0,
    n_jobs: int = 1,
) -> np.ndarray:
    """
    Extract locus sub-regions from each tile in a tiled chromosome simulation,
    compute contact maps per tile, and average across all tiles and frames.

    This implements the "tiling trick" from the Hansen lab: a single simulation
    of a large chromosome contains multiple identical copies of the locus
    region. By extracting and averaging contacts from all copies, we get N×
    more data per simulation frame (N = n_tiles from TILING config).

    Memory-efficient: tiles are yielded lazily one at a time from the
    conformations iterable — the full tile list is never materialised in RAM.
    With 4000 frames × 28 tiles this saves ~5.4 GB compared to building a list.

    Parameters
    ----------
    conformations : list or iterable of np.ndarray
        Each element is (chrom_size, 3) coordinates from a tiled simulation.
    contact_radius : float
        Distance threshold for contact detection.
    n_jobs : int
        Number of parallel workers.

    Returns
    -------
    contact_map : np.ndarray of shape (N_MONOMERS, N_MONOMERS)
        Averaged contact map for the locus region.
    """
    n_tiles  = TILING["n_tiles"]
    tile_size = TILING["tile_size"]
    pad      = TILING["padding"]
    N_locus  = N_MONOMERS  # 2000
    tile_offsets = [t * tile_size + pad for t in range(n_tiles)]

    # ------------------------------------------------------------------
    # Preferred path: file-parallel (Hansen-lab pattern).
    # When the input carries a list of polychrom blocks_*.h5 paths, hand
    # the file list to workers and let them read their own data. Parent
    # heap stays ~0 regardless of trajectory length — no frame ever lives
    # in the parent, so fork() CoW inheritance can't amplify.
    # ------------------------------------------------------------------
    block_files = getattr(conformations, "block_files", None)
    if block_files and n_jobs > 1:
        # We still need chrom_size to decide tiled vs non-tiled. Peek at
        # the first block's first frame WITHOUT consuming the stream, so
        # the worker pool sees the full file list.
        import h5py as _h5py
        with _h5py.File(block_files[0], "r") as _hf:
            _keys = sorted(_hf.keys(), key=lambda x: int(x) if x.isdigit() else 0)
            _first_key = next(k for k in _keys
                              if isinstance(_hf[k], _h5py.Group) and "pos" in _hf[k])
            chrom_size = int(_hf[_first_key]["pos"].shape[0])
        if chrom_size == N_locus:
            # Single-locus: one "tile" at offset 0.
            return _extract_tiles_file_parallel(
                block_files, contact_radius, n_jobs,
                tile_offsets=[0], N_locus=N_locus,
            )
        if chrom_size == TILING["chrom_size"]:
            return _extract_tiles_file_parallel(
                block_files, contact_radius, n_jobs,
                tile_offsets=tile_offsets, N_locus=N_locus,
            )
        logger.warning(
            f"Unexpected conformation size {chrom_size}; block_files present "
            f"but neither N_locus={N_locus} nor TILING['chrom_size']="
            f"{TILING['chrom_size']}. Falling back to streaming path."
        )

    # ------------------------------------------------------------------
    # Streaming fallback (used when block_files is unknown, e.g. in-memory
    # lists handed in by tests).
    # ------------------------------------------------------------------
    conf_iter = iter(conformations)
    try:
        first_frame = next(conf_iter)
    except StopIteration:
        raise RuntimeError("No conformations to process.")

    chrom_size = first_frame.shape[0]

    def _full_iter():
        """Re-yield the first frame then the rest."""
        yield first_frame
        yield from conf_iter

    if chrom_size == N_locus:
        logger.info("Conformations are single-locus (not tiled), computing directly...")
        cmap, n = _accumulate_contacts(N_locus, contact_radius, _full_iter(), n_jobs)
        return cmap / n

    if chrom_size != TILING["chrom_size"]:
        logger.warning(f"Unexpected conformation size {chrom_size}, expected "
                       f"{TILING['chrom_size']} (tiled) or {N_locus} (single). "
                       f"Computing on full size.")
        cmap, n = _accumulate_contacts(chrom_size, contact_radius, _full_iter(), n_jobs)
        return cmap / n

    n_frames = len(conformations) if hasattr(conformations, '__len__') else None
    if n_frames is not None:
        total_tiles = n_frames * n_tiles
        logger.info(f"Extracting {n_tiles} tiles/frame from {n_frames} frames "
                    f"= {total_tiles:,} sub-regions (streamed, no tile list built)")
    else:
        logger.info(f"Extracting {n_tiles} tiles/frame (streaming)...")

    def _tile_iter():
        """Generator: yield each locus sub-region from every frame in order."""
        for frame_coords in _full_iter():
            for t in range(n_tiles):
                start = t * tile_size + pad
                yield frame_coords[start:start + N_locus]

    total_tiles_hint = total_tiles if n_frames is not None else None
    cmap, n = _accumulate_contacts(N_locus, contact_radius, _tile_iter(), n_jobs,
                                   n_frames=total_tiles_hint)
    return cmap / n


class ConformationStream:
    """Lazy iterable over polymer conformations stored on disk.

    Memory footprint is ~one frame (a few MB) regardless of trajectory
    length. Supersedes the old list-materialising loader which reached
    ~170 GB for 100k-frame tiled (70k-monomer) trajectories and caused
    OOM during contact-map tile extraction.

    Supported operations:
        len(stream)       -> total frame count (known upfront from metadata)
        iter(stream)      -> fresh generator each call; streams frames from
                             disk one at a time. Safe to iterate multiple
                             times (e.g. per-module passes for contact_maps,
                             polymer_dynamics, msd_two_point), each pass
                             re-reads the file.
        stream[0]         -> cached first frame. The consumers in
                             polymer_dynamics / msd_two_point use this only
                             to peek .shape[0] (chromosome size).

    Arbitrary indexing (stream[17]) is intentionally unsupported; it would
    silently re-enable the old load-everything pattern. Iterate or
    materialise a slice explicitly if you need random access.
    """
    __slots__ = ("_n", "_iter_factory", "_source_label", "_head", "block_files")

    def __init__(self, *, n: int, iter_factory, source_label: str,
                 block_files: Optional[list] = None):
        self._n = n
        self._iter_factory = iter_factory
        self._source_label = source_label
        self._head = None
        # If non-None, the concrete list of polychrom blocks_*.h5 paths backing
        # this stream. Lets `extract_tiles_and_average` hand filenames (not
        # frames) to workers so the parent heap never touches trajectory data.
        # See the Hansen-lab pattern in polychrom/contactmaps.py:averageContacts.
        self.block_files = block_files

    def __len__(self):
        return self._n

    def __iter__(self):
        return self._iter_factory()

    def __getitem__(self, idx):
        if idx == 0 or idx == -self._n:
            if self._head is None:
                self._head = next(iter(self))
            return self._head
        raise IndexError(
            f"ConformationStream[{idx}] not supported — only [0] is cached. "
            f"Iterate with `for f in stream`, or materialise with list(stream).")

    def __repr__(self):
        return f"<ConformationStream n={self._n} source={self._source_label!r}>"

    @classmethod
    def chain(cls, *streams) -> "ConformationStream":
        """Chain multiple streams end-to-end into a single re-iterable stream.

        Used by Phase 2 pooling to pass conformations from N replicates to
        `extract_tiles_and_average` and `_run_msd_and_dynamics` without
        materialising them into a list (which would OOM at ~500 GB for
        3 × 100k-frame trajectories).

        Each call to iter() re-reads from disk for every underlying stream.
        If every input stream carries a non-None ``block_files`` list, the
        chained stream merges them so the file-parallel contact-map path
        works on pooled conditions too.
        """
        n = sum(len(s) for s in streams)

        def _iter():
            for s in streams:
                yield from s

        merged_blocks: Optional[list] = None
        per_stream_blocks = [getattr(s, "block_files", None) for s in streams]
        if all(bf is not None for bf in per_stream_blocks):
            merged_blocks = []
            for bf in per_stream_blocks:
                merged_blocks.extend(bf)

        return cls(n=n, iter_factory=_iter,
                   source_label=f"chain<{len(streams)}>",
                   block_files=merged_blocks)


def load_conformations_h5(sim_dir: str):
    """
    Open conformations from a simulation output directory as a lazy stream.

    Returns a `ConformationStream` that reads one frame at a time from disk
    on each iteration. The in-memory footprint is a few MB (one frame),
    not the full trajectory.

    Supports the same on-disk formats as before:
      1. Polychrom HDF5Reporter format (blocks_*.h5) — preferred
      2. Raw conformations.h5 (frame_0, frame_1, ...) — legacy/OpenMM fallback
      3. Returns None if only LEF-only contact maps are present
    """
    import h5py
    import glob

    # --- Polychrom HDF5Reporter format (blocks_*.h5) ---
    # We intentionally bypass polychrom.hdf5_format.load_URI here: that helper
    # opens and closes each HDF5 file per frame (100k opens for a 100k-frame
    # trajectory), which is both slow and — more importantly — appeared to
    # prevent the HDF5 library cache from being released between calls,
    # producing apparent "memory leaks" during streaming tile extraction.
    # By iterating each block file exactly once (open → read all its frames →
    # close), we get deterministic resource release and far less overhead.
    block_files = sorted(glob.glob(os.path.join(sim_dir, "blocks_*.h5")))
    if block_files:
        entries = []
        for bf in block_files:
            with h5py.File(bf, "r") as hf:
                for key in sorted(hf.keys(), key=lambda x: int(x) if x.isdigit() else 0):
                    if isinstance(hf[key], h5py.Group) and "pos" in hf[key]:
                        entries.append((bf, key))
        if entries:
            def _iter():
                cur_bf = None
                cur_hf = None
                try:
                    for bf, key in entries:
                        if bf != cur_bf:
                            if cur_hf is not None:
                                cur_hf.close()
                            cur_hf = h5py.File(bf, "r")
                            cur_bf = bf
                        # np.asarray(..., dtype=...) forces a contiguous copy
                        # into a fresh numpy buffer, decoupling the yielded
                        # frame from the HDF5 dataset object.
                        yield np.asarray(cur_hf[key]["pos"], dtype=np.float64)
                finally:
                    if cur_hf is not None:
                        cur_hf.close()
            logger.info(f"Opened {len(entries)} conformations (streaming, direct h5py) "
                        f"from {len(block_files)} block files")
            return ConformationStream(n=len(entries), iter_factory=_iter,
                                      source_label="direct_h5py_blocks",
                                      block_files=list(block_files))

    # --- Try raw conformations.h5 (OpenMM standalone format) ---
    h5_path = os.path.join(sim_dir, "conformations.h5")
    if os.path.exists(h5_path):
        with h5py.File(h5_path, "r") as hf:
            n = int(hf.attrs.get("n_frames",
                                 len([k for k in hf.keys() if k.startswith("frame_")])))
        def _iter():
            with h5py.File(h5_path, "r") as hf:
                for i in range(n):
                    yield hf[f"frame_{i}"][:]
        logger.info(f"Opened {n} conformations (streaming) from {h5_path}")
        return ConformationStream(n=n, iter_factory=_iter,
                                  source_label="legacy_h5")

    # --- LEF-only fallback ---
    lef_map_path = os.path.join(sim_dir, "lef_contact_map.npy")
    if os.path.exists(lef_map_path):
        logger.info("No 3D conformations found. Using LEF bridging contact map.")
        return None

    raise FileNotFoundError(f"No conformation data found in {sim_dir}")


def load_lef_contact_map(sim_dir: str) -> Optional[np.ndarray]:
    """Load pre-computed LEF bridging contact map."""
    path = os.path.join(sim_dir, "lef_contact_map.npy")
    if os.path.exists(path):
        return np.load(path)
    return None


# =============================================================================
# CONTACT PROBABILITY P(s)
# =============================================================================

def compute_ps_curve(contact_map: np.ndarray, resolution: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute P(s) — contact probability as a function of genomic distance.

    Vectorised via np.diagonal: each diagonal of the contact map corresponds
    to a fixed genomic separation, so we just average each diagonal.
    This replaces an O(N²) Python double loop with N fast C-level mean calls
    — ~40× faster for N=2000.

    Parameters
    ----------
    contact_map : np.ndarray
        Square contact frequency matrix.
    resolution : int
        Bp per bin.

    Returns
    -------
    distances : np.ndarray
        Genomic distances in bp.
    ps : np.ndarray
        Average contact probability at each distance.
    """
    N = contact_map.shape[0]
    # np.diagonal(k) returns a read-only view of the k-th diagonal — no copy.
    ps = np.array([np.mean(np.diagonal(contact_map, offset=d)) for d in range(N)])
    distances = np.arange(N) * resolution
    return distances, ps


# =============================================================================
# P(s) POWER LAW FITTING (Hansen lab style)
# =============================================================================

def fit_ps_powerlaw(
    distances: np.ndarray,
    ps: np.ndarray,
    lower_bound_bp: int = 2000,
    upper_bound_bp: int = 30000,
) -> dict:
    """
    Fit P(s) with a piecewise power law in log-log space.

    Mirny-lab / Gassler-style log-log fit (NOT a direct port of any
    Hansen-lab routine — their ``looptools.calculate_and_save_avg_Ps_curve``
    in github.com/ahansenlab/AbsQuant_analysis_code wraps
    ``cooltools.expected_cis`` and does not extract a power-law exponent):
      - Region 1 (s < lower_bound): 6th-degree polynomial in log-log space
        to capture the non-scaling short-range behaviour. Note: at the
        default 1 kb resolution this region has too few bins (≤ 2) for
        the deg-6 polyfit's ``>7`` requirement and is silently skipped.
      - Region 2 (lower_bound <= s < upper_bound): 1st-degree polynomial
        in log-log (= true power law, P(s) ~ s^(-alpha)).

    The power law exponent alpha is the key biophysical parameter: it
    reflects the polymer state (alpha ~ 1.0-1.1 for fractal globule,
    alpha ~ 1.5 for equilibrium globule, alpha ~ 1.3 for typical
    mammalian chromatin with loop extrusion).

    Parameters
    ----------
    distances : np.ndarray
        Genomic distances in bp (from compute_ps_curve).
    ps : np.ndarray
        Contact probabilities at each distance.
    lower_bound_bp : int
        Transition from polynomial to power law (default 2 kb).
    upper_bound_bp : int
        Upper limit of the power-law fit (default 30 kb).

    Returns
    -------
    fit : dict
        'exponent'       : float — power law exponent alpha (P(s) ~ s^-alpha)
        'log_intercept'  : float — intercept of the log-log linear fit
        'coeffs_region1' : np.ndarray — polynomial coefficients (short range)
        'coeffs_region2' : np.ndarray — [slope, intercept] of log-log linear fit
        'lower_bound_bp' : int
        'upper_bound_bp' : int
    """
    # Filter out zeros and NaN for log-space fitting
    valid = (distances > 0) & (ps > 0) & np.isfinite(ps)
    d_valid = distances[valid]
    p_valid = ps[valid]

    log_d = np.log10(d_valid)
    log_p = np.log10(p_valid)

    fit = {
        "lower_bound_bp": lower_bound_bp,
        "upper_bound_bp": upper_bound_bp,
    }

    # Region 1: short range — 6th-degree polynomial in log-log
    mask1 = d_valid < lower_bound_bp
    if mask1.sum() > 7:
        fit["coeffs_region1"] = np.polyfit(log_d[mask1], log_p[mask1], deg=6).tolist()
    else:
        fit["coeffs_region1"] = None

    # Region 2: power-law range — linear fit in log-log
    mask2 = (d_valid >= lower_bound_bp) & (d_valid < upper_bound_bp)
    if mask2.sum() > 2:
        coeffs2 = np.polyfit(log_d[mask2], log_p[mask2], deg=1)
        fit["coeffs_region2"] = coeffs2.tolist()
        fit["exponent"] = float(-coeffs2[0])       # alpha in P(s) ~ s^-alpha
        fit["log_intercept"] = float(coeffs2[1])
    else:
        fit["coeffs_region2"] = None
        fit["exponent"] = None
        fit["log_intercept"] = None

    return fit


def extract_ps_metrics(
    distances: np.ndarray,
    ps: np.ndarray,
    resolution: int = 1000,
) -> dict:
    """
    Extract key metrics from a P(s) curve.

    Returns contact probabilities at standard genomic distances, the
    power law exponent, and the log-log derivative (slope) at several
    distances.

    Parameters
    ----------
    distances : np.ndarray
        Genomic distances in bp.
    ps : np.ndarray
        Contact probabilities.
    resolution : int
        Bp per bin.

    Returns
    -------
    metrics : dict
        Contact probabilities at key distances, power-law exponent,
        and log-log derivative at selected points.
    """
    metrics = {}

    # ── Contact probability at standard distances ────────────────────
    for label, target_bp in [("10kb", 10_000), ("50kb", 50_000),
                             ("100kb", 100_000), ("200kb", 200_000),
                             ("500kb", 500_000), ("1Mb", 1_000_000)]:
        idx = int(round(target_bp / resolution))
        if 0 <= idx < len(ps):
            metrics[f"P(s={label})"] = float(ps[idx])

    # ── Power law fit ────────────────────────────────────────────────
    fit = fit_ps_powerlaw(distances, ps)
    if fit["exponent"] is not None:
        metrics["ps_exponent_2_30kb"] = fit["exponent"]
        metrics["ps_log_intercept"] = fit["log_intercept"]

    # ── Local log-log slope at several distances ─────────────────────
    # dlog(P)/dlog(s) computed as finite difference on the log-log curve.
    # This is the "running exponent" — it reveals where the power law
    # holds (constant slope) and where it breaks down.
    valid = (distances > 0) & (ps > 0) & np.isfinite(ps)
    d_v = distances[valid]
    p_v = ps[valid]
    if len(d_v) > 5:
        log_d = np.log10(d_v)
        log_p = np.log10(p_v)
        # Central difference for derivative
        dlogp_dlogs = np.gradient(log_p, log_d)
        for label, target_bp in [("5kb", 5_000), ("10kb", 10_000),
                                 ("50kb", 50_000), ("100kb", 100_000),
                                 ("500kb", 500_000)]:
            idx = np.argmin(np.abs(d_v - target_bp))
            if 2 <= idx < len(dlogp_dlogs) - 2:
                metrics[f"dlogP_dlogs_at_{label}"] = float(dlogp_dlogs[idx])

    return metrics


# =============================================================================
# INSULATION SCORE
# =============================================================================

def compute_insulation_score(
    contact_map: np.ndarray,
    window: int = 50,
) -> np.ndarray:
    """
    Compute diamond insulation score (Crane et al. 2015).

    For each position i, computes the mean contact frequency in the
    off-diagonal block [i-w:i, i:i+w] relative to the flanking on-diagonal
    blocks [i-w:i, i-w:i] and [i:i+w, i:i+w].

    Parameters
    ----------
    contact_map : np.ndarray
        Square contact matrix.
    window : int
        Window size in bins for the diamond.

    Returns
    -------
    insulation : np.ndarray
        Log2 insulation score at each position (NaN at boundaries).
    """
    N = contact_map.shape[0]
    insulation = np.full(N, np.nan)

    for i in range(window, N - window):
        # Sum contacts in the off-diagonal block crossing position i
        diamond      = contact_map[i - window:i, i:i + window]
        diamond_mean = np.nanmean(diamond)

        if diamond_mean > 0:
            # Normalize by mean of flanking on-diagonal blocks
            left_diamond  = contact_map[i - window:i, i - window:i]
            right_diamond = contact_map[i:i + window, i:i + window]
            flank_mean    = (np.nanmean(left_diamond) + np.nanmean(right_diamond)) / 2

            if flank_mean > 0:
                insulation[i] = np.log2(diamond_mean / flank_mean)

    return insulation


# =============================================================================
# COMPARISON METRICS
# =============================================================================

def compare_contact_maps(
    sim_map: np.ndarray,
    exp_map: np.ndarray,
    max_diag: Optional[int] = None,
) -> dict:
    """
    Compare simulated and experimental contact maps using multiple metrics.

    Parameters
    ----------
    sim_map : np.ndarray
        Simulated contact matrix.
    exp_map : np.ndarray
        Experimental contact matrix.
    max_diag : int, optional
        Maximum diagonal distance to consider.

    Returns
    -------
    metrics : dict
        Dictionary of comparison metrics.
    """
    N = min(sim_map.shape[0], exp_map.shape[0])
    sim = sim_map[:N, :N].copy()
    exp = exp_map[:N, :N].copy()

    if max_diag is not None:
        # Vectorised masking — replaces O(N²) Python loop
        row_idx, col_idx = np.ogrid[:N, :N]
        off_diag_mask = np.abs(row_idx - col_idx) > max_diag
        sim[off_diag_mask] = 0
        exp[off_diag_mask] = 0

    # Normalize both maps to same scale
    sim_norm = sim / (sim.max() + 1e-10)
    exp_norm = exp / (exp.max() + 1e-10)

    # 1. Pearson correlation (per-diagonal)
    # np.diagonal(offset=d) returns a zero-copy view of diagonal d — no list comp needed.
    diag_corrs = []
    for d in range(1, min(N, max_diag or N)):
        sim_diag = np.diagonal(sim_norm, offset=d)
        exp_diag = np.diagonal(exp_norm, offset=d)
        if np.std(sim_diag) > 0 and np.std(exp_diag) > 0:
            corr = np.corrcoef(sim_diag, exp_diag)[0, 1]
            diag_corrs.append(corr)

    # 2. Overall Pearson correlation (upper triangle)
    upper_idx = np.triu_indices(N, k=1)
    sim_upper = sim_norm[upper_idx]
    exp_upper = exp_norm[upper_idx]
    mask = (sim_upper > 0) | (exp_upper > 0)
    if mask.sum() > 0:
        overall_corr = np.corrcoef(sim_upper[mask], exp_upper[mask])[0, 1]
    else:
        overall_corr = 0.0

    # 3. Stratum-adjusted correlation (SCC-like)
    # Simple version: weighted average of per-diagonal correlations
    weights = np.array([N - d for d in range(1, len(diag_corrs) + 1)], dtype=float)
    weights /= weights.sum()
    scc = np.sum(np.array(diag_corrs[:len(weights)]) * weights) if diag_corrs else 0.0

    # 4. P(s) curve comparison
    _, ps_sim = compute_ps_curve(sim)
    _, ps_exp = compute_ps_curve(exp)
    min_len = min(len(ps_sim), len(ps_exp))
    ps_mask = (ps_sim[:min_len] > 0) & (ps_exp[:min_len] > 0)
    if ps_mask.sum() > 2:
        ps_corr = np.corrcoef(
            np.log10(ps_sim[:min_len][ps_mask] + 1e-10),
            np.log10(ps_exp[:min_len][ps_mask] + 1e-10)
        )[0, 1]
    else:
        ps_corr = 0.0

    # 5. Insulation score correlation
    ins_sim = compute_insulation_score(sim)
    ins_exp = compute_insulation_score(exp)
    valid = ~(np.isnan(ins_sim) | np.isnan(ins_exp))
    if valid.sum() > 10:
        ins_corr = np.corrcoef(ins_sim[valid], ins_exp[valid])[0, 1]
    else:
        ins_corr = 0.0

    # 6. P(s) power law exponents for both sim and exp
    sim_ps_metrics = extract_ps_metrics(np.arange(N) * 1000, ps_sim[:N])
    exp_ps_metrics = extract_ps_metrics(np.arange(N) * 1000, ps_exp[:N])

    metrics = {
        "overall_pearson": float(overall_corr),
        "mean_diag_pearson": float(np.nanmean(diag_corrs)) if diag_corrs else 0.0,
        "stratum_adjusted_corr": float(scc),
        "ps_curve_corr": float(ps_corr),
        "insulation_corr": float(ins_corr),
        "n_diagonals_compared": len(diag_corrs),
        # P(s) power law exponents
        "sim_ps_exponent_2_30kb": sim_ps_metrics.get("ps_exponent_2_30kb"),
        "exp_ps_exponent_2_30kb": exp_ps_metrics.get("ps_exponent_2_30kb"),
    }

    # Add per-distance contact probabilities for sim and exp
    for key, val in sim_ps_metrics.items():
        if key.startswith("P(s=") or key.startswith("dlogP_"):
            metrics[f"sim_{key}"] = val
    for key, val in exp_ps_metrics.items():
        if key.startswith("P(s=") or key.startswith("dlogP_"):
            metrics[f"exp_{key}"] = val

    return metrics


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Compute and compare contact maps")
    parser.add_argument("--sim-dir", type=str, required=True,
                        help="Simulation results directory")
    parser.add_argument("--exp-hic", type=str, default=None,
                        help="Path to experimental Hi-C .npy matrix")
    parser.add_argument("--contact-radius", type=float, default=SIM_RUN["contact_radius"],
                        help="Contact detection radius")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory for analysis results")
    parser.add_argument("--n-jobs", type=int, default=4,
                        help="Number of parallel workers for contact map computation")

    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(args.sim_dir, "analysis")
    os.makedirs(args.output, exist_ok=True)

    # --- Load or compute simulated contact map ---
    conformations = None
    try:
        conformations = load_conformations_h5(args.sim_dir)
    except FileNotFoundError:
        pass

    if conformations is not None:
        n = len(conformations) if hasattr(conformations, '__len__') else '?'
        logger.info(f"Computing contact map from {n} conformations...")
        # Use tile extraction if conformations are from a tiled simulation
        sim_map = extract_tiles_and_average(
            conformations, args.contact_radius, n_jobs=args.n_jobs)
    else:
        sim_map = load_lef_contact_map(args.sim_dir)
        if sim_map is None:
            raise RuntimeError(f"No simulation data found in {args.sim_dir}")
        logger.info("Using LEF bridging contact map")

    # Save simulated contact map
    np.save(os.path.join(args.output, "sim_contact_map.npy"), sim_map)

    # P(s) curve
    distances, ps = compute_ps_curve(sim_map, RESOLUTION)
    np.savez(os.path.join(args.output, "sim_ps_curve.npz"),
             distances=distances, ps=ps)

    # Insulation score
    insulation = compute_insulation_score(sim_map)
    np.save(os.path.join(args.output, "sim_insulation.npy"), insulation)

    logger.info(f"Simulated contact map saved (shape: {sim_map.shape})")

    # --- Compare with experimental Hi-C ---
    if args.exp_hic and os.path.exists(args.exp_hic):
        logger.info("Comparing with experimental Hi-C...")
        exp_map = np.load(args.exp_hic)
        metrics = compare_contact_maps(sim_map, exp_map)

        with open(os.path.join(args.output, "comparison_metrics.json"), "w") as f:
            json.dump(metrics, f, indent=2)

        logger.info("Comparison metrics:")
        for k, v in metrics.items():
            logger.info(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
    else:
        logger.info("No experimental Hi-C provided. Skipping comparison.")
        logger.info("To compare, provide --exp-hic path to .npy contact matrix")

    logger.info(f"Analysis results saved to: {args.output}")


if __name__ == "__main__":
    main()
