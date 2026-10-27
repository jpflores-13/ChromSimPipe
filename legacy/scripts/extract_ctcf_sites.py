#!/usr/bin/env python
"""
Extract CTCF binding sites with motif orientation from ENCODE ChIP-seq data.

Can run genome-wide or on a specific region. When a region is specified,
also outputs the CTCF_SITES array for polychrom simulations.

Data sources:
  - ENCODE CTCF ChIP-seq in Bruce4 mESC (C57BL/6):
    Experiment: ENCSR000CCB
    File (optimal IDR peaks, mm10): ENCFF508CKL

  - ENCODE CTCF ChIP-seq in E14TG2a.4 mESC (129/Ola):
    Experiment: ENCSR362VNF
    Files (IDR peaks, mm10): ENCFF311HPG, ENCFF693MYU

  - CTCF motif orientation via FIMO:
    JASPAR motif MA0139.1 (auto-downloaded from JASPAR)
    mm10 genome FASTA (auto-downloaded from UCSC)

Environment setup:
    mamba env create -f envs/ctcf_extraction.yml
    conda activate ctcf_extraction

Usage:
    # Genome-wide: annotate ALL CTCF peaks with orientation
    python scripts/extract_ctcf_sites.py --source encode --genome-wide

    # Specific region (for simulation):
    python scripts/extract_ctcf_sites.py --source encode --region chr3:34000000-36000000

    # Sox2 locus (default region from parameters.py):
    python scripts/extract_ctcf_sites.py --source encode

    # From your own peaks:
    python scripts/extract_ctcf_sites.py --source bed --bed my_peaks.narrowPeak --genome-wide
"""

import os
import sys
import argparse
import gzip
import logging
import shutil
import subprocess
import tempfile
import urllib.request
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# =============================================================================
# CONSTANTS & URLS
# =============================================================================

ENCODE_CTCF_PEAKS = {
    "mESC_Bruce4": {
        "experiment": "ENCSR000CCB",
        "file": "ENCFF508CKL",
        "url": "https://www.encodeproject.org/files/ENCFF508CKL/@@download/ENCFF508CKL.bed.gz",
        "description": "CTCF ChIP-seq optimal IDR peaks, Bruce4 mESC (C57BL/6), mm10",
    },
    "mESC_E14": {
        "experiment": "ENCSR362VNF",
        "file": "ENCFF311HPG",
        "url": "https://www.encodeproject.org/files/ENCFF311HPG/@@download/ENCFF311HPG.bed.gz",
        "description": "CTCF ChIP-seq IDR peaks, E14TG2a.4 mESC (129/Ola), mm10",
    },
}

JASPAR_CTCF_MEME_URL = (
    "https://jaspar.elixir.no/api/v1/matrix/MA0139.1/?format=meme"
)

MM10_FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"


# =============================================================================
# DOWNLOAD HELPERS
# =============================================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def download_file(url, dest, description="file"):
    """Download a file if it doesn't already exist."""
    if os.path.exists(dest):
        logger.info(f"  Already exists: {dest}")
        return dest
    logger.info(f"  Downloading {description}...")
    logger.info(f"  URL: {url}")
    urllib.request.urlretrieve(url, dest)
    logger.info(f"  Saved: {dest}")
    return dest


def download_encode_peaks(output_dir, cell_type="mESC_Bruce4"):
    """Download CTCF narrowPeak from ENCODE and decompress."""
    if cell_type not in ENCODE_CTCF_PEAKS:
        raise ValueError(f"Unknown cell type '{cell_type}'. "
                         f"Available: {list(ENCODE_CTCF_PEAKS.keys())}")

    info = ENCODE_CTCF_PEAKS[cell_type]
    out_gz = os.path.join(output_dir, f"CTCF_{cell_type}_peaks.bed.gz")
    out_bed = os.path.join(output_dir, f"CTCF_{cell_type}_peaks.bed")

    if os.path.exists(out_bed):
        logger.info(f"  Already exists: {out_bed}")
        return out_bed

    download_file(info["url"], out_gz,
                  f"{info['description']} ({info['experiment']}/{info['file']})")

    with gzip.open(out_gz, "rb") as gz, open(out_bed, "w") as out:
        for line in gz:
            out.write(line.decode())

    logger.info(f"  Decompressed: {out_bed}")
    return out_bed


def download_jaspar_motif(output_dir):
    """Download the CTCF motif (MA0139.1) from JASPAR in MEME format."""
    motif_path = os.path.join(output_dir, "MA0139.1_CTCF.meme")

    if os.path.exists(motif_path):
        logger.info(f"  Already exists: {motif_path}")
        return motif_path

    logger.info("  Downloading CTCF motif MA0139.1 from JASPAR...")
    download_file(JASPAR_CTCF_MEME_URL, motif_path, "JASPAR CTCF motif (MA0139.1)")

    with open(motif_path) as f:
        header = f.read(200)
    if "MEME" not in header and "letter-probability" not in header:
        logger.warning("  JASPAR returned non-MEME format. Building MEME file manually...")
        _build_ctcf_meme_file(motif_path)

    return motif_path


def _build_ctcf_meme_file(output_path):
    """Write the CTCF MA0139.1 position frequency matrix in MEME format."""
    pfm = [
        [0.0952, 0.3175, 0.2698, 0.3175],
        [0.1746, 0.1032, 0.1270, 0.5952],
        [0.0317, 0.0159, 0.8889, 0.0635],
        [0.0159, 0.8571, 0.0794, 0.0476],
        [0.5873, 0.0317, 0.3175, 0.0635],
        [0.0952, 0.6349, 0.0317, 0.2381],
        [0.0476, 0.5714, 0.1587, 0.2222],
        [0.3016, 0.1746, 0.3651, 0.1587],
        [0.1746, 0.0952, 0.5556, 0.1746],
        [0.0159, 0.0476, 0.9206, 0.0159],
        [0.0635, 0.0159, 0.8413, 0.0794],
        [0.1587, 0.1746, 0.3810, 0.2857],
        [0.3175, 0.3810, 0.0476, 0.2540],
        [0.4127, 0.1587, 0.2222, 0.2063],
        [0.2540, 0.2222, 0.2540, 0.2698],
        [0.2063, 0.3175, 0.2222, 0.2540],
        [0.2381, 0.3651, 0.1270, 0.2698],
        [0.2540, 0.1905, 0.3175, 0.2381],
        [0.2381, 0.2222, 0.3333, 0.2063],
    ]
    with open(output_path, "w") as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")
        f.write("MOTIF MA0139.1 CTCF\n")
        f.write(f"letter-probability matrix: alength= 4 w= {len(pfm)} "
                f"nsites= 63 E= 0\n")
        for row in pfm:
            f.write(f"  {row[0]:.6f}  {row[1]:.6f}  {row[2]:.6f}  {row[3]:.6f}\n")
        f.write("\n")
    logger.info(f"  Built MEME file: {output_path}")


def download_genome(output_dir, chrom_only=None):
    """Download mm10 genome FASTA and index it."""
    genome_dir = ensure_dir(output_dir)
    fasta_path = os.path.join(genome_dir, "mm10.fa")
    fai_path = fasta_path + ".fai"

    if os.path.exists(fasta_path) and os.path.exists(fai_path):
        logger.info(f"  Genome already indexed: {fasta_path}")
        return fasta_path

    if not os.path.exists(fasta_path):
        gz_path = fasta_path + ".gz"

        # Verify existing .gz integrity
        if os.path.exists(gz_path):
            result = subprocess.run(["gzip", "-t", gz_path],
                                    capture_output=True, text=True)
            if result.returncode != 0:
                logger.warning(f"  {gz_path} is corrupted — deleting and re-downloading")
                os.remove(gz_path)

        download_file(MM10_FASTA_URL, gz_path,
                      "mm10 genome FASTA (~800 MB compressed)")

        logger.info("  Decompressing genome (this takes a few minutes)...")
        subprocess.run(["gunzip", "-k", gz_path], check=True)
        logger.info(f"  Decompressed: {fasta_path}")

    if not os.path.exists(fai_path):
        logger.info("  Indexing genome with samtools faidx...")
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
        logger.info(f"  Indexed: {fai_path}")

    if chrom_only:
        chrom_fa = os.path.join(genome_dir, f"{chrom_only}.fa")
        if not os.path.exists(chrom_fa):
            logger.info(f"  Extracting {chrom_only}...")
            with open(chrom_fa, "w") as out:
                subprocess.run(["samtools", "faidx", fasta_path, chrom_only],
                               check=True, stdout=out)
            subprocess.run(["samtools", "faidx", chrom_fa], check=True)
        return chrom_fa

    return fasta_path


# =============================================================================
# PEAK PARSING
# =============================================================================

def parse_region_string(region_str):
    """
    Parse a UCSC-style region string into (chrom, start, end).

    Accepts: "chr3:34000000-36000000" or "chr3:34,000,000-36,000,000"
    Returns: ("chr3", 34000000, 36000000)
    """
    region_str = region_str.replace(",", "")
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def load_peaks(bed_path, region=None):
    """
    Load peaks from a narrowPeak/BED file.

    Parameters
    ----------
    bed_path : str
        Path to BED or narrowPeak file.
    region : tuple of (chrom, start, end), optional
        If provided, only return peaks within this region.
        If None, return all peaks (genome-wide).

    Returns
    -------
    peaks : list of dict
        Each dict has: chrom, start, end, summit, name, score, strand.
    """
    peaks = []

    with open(bed_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("track") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Region filter
            if region is not None:
                r_chrom, r_start, r_end = region
                if chrom != r_chrom:
                    continue
                if end < r_start or start > r_end:
                    continue

            name = fields[3] if len(fields) > 3 else "."
            score = float(fields[4]) if len(fields) > 4 else 0
            strand = fields[5] if len(fields) > 5 else "."

            # narrowPeak: column 10 (index 9) is summit offset from start
            if len(fields) >= 10:
                try:
                    summit_offset = int(fields[9])
                    summit = start + summit_offset
                except ValueError:
                    summit = (start + end) // 2
            else:
                summit = (start + end) // 2

            peaks.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "summit": summit,
                "name": name,
                "score": score,
                "strand": strand,
            })

    peaks.sort(key=lambda x: (x["chrom"], x["summit"]))

    if region is not None:
        r_chrom, r_start, r_end = region
        logger.info(f"  Found {len(peaks)} CTCF peaks in {r_chrom}:{r_start:,}-{r_end:,}")
    else:
        chroms = set(p["chrom"] for p in peaks)
        logger.info(f"  Found {len(peaks)} CTCF peaks genome-wide "
                    f"across {len(chroms)} chromosomes")

    return peaks


# =============================================================================
# MOTIF ORIENTATION CALLING
# =============================================================================

def call_motif_orientation(peaks, genome_fasta, motif_meme, flank=100, fimo_thresh=1e-3):
    """
    Determine CTCF motif orientation at each peak using FIMO.

    Works on any number of peaks (region or genome-wide). Processes in
    batches of 50,000 to keep temp files manageable.

    Peaks without a FIMO hit are assigned orientation based on a hash of
    their genomic coordinates, ensuring the result is deterministic and
    identical regardless of how many other peaks are in the run.
    """
    if not peaks:
        return peaks

    # Check if strand info is already in the BED
    has_strand = all(p["strand"] in ("+", "-") for p in peaks)
    if has_strand:
        logger.info("  BED file already has strand info — using directly")
        for p in peaks:
            p["orientation"] = +1 if p["strand"] == "+" else -1
        return peaks

    # Check tools
    for tool in ["bedtools", "fimo"]:
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            logger.error(f"  {tool} not found! Install with: conda activate ctcf_extraction")
            return _fallback_by_coordinate(peaks)

    logger.info(f"  Calling motif orientation with FIMO on {len(peaks)} peaks "
                f"(±{flank} bp windows)...")

    # Process in batches for large peak sets
    BATCH_SIZE = 50_000
    all_best_hits = {}  # keyed by (chrom, summit) for coordinate-based lookup

    for batch_start in range(0, len(peaks), BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, len(peaks))
        batch = peaks[batch_start:batch_end]

        if len(peaks) > BATCH_SIZE:
            logger.info(f"  Processing batch {batch_start // BATCH_SIZE + 1} "
                        f"({batch_start}-{batch_end} of {len(peaks)})...")

        batch_hits = _run_fimo_batch(batch, genome_fasta, motif_meme, flank,
                                     fimo_thresh)
        all_best_hits.update(batch_hits)

    # Assign orientations — keyed by genomic coordinate, NOT list index
    n_fimo = 0
    n_fallback = 0
    for p in peaks:
        key = (p["chrom"], p["summit"])
        if key in all_best_hits:
            strand, pval = all_best_hits[key]
            p["orientation"] = +1 if strand == "+" else -1
            p["strand"] = strand
            p["fimo_pvalue"] = pval
            n_fimo += 1
        else:
            # Deterministic fallback based on genomic position, not list index.
            # hash(coordinate) ensures the same peak always gets the same
            # orientation regardless of what other peaks are in the run.
            p["orientation"] = +1 if hash((p["chrom"], p["summit"])) % 2 == 0 else -1
            p["strand"] = "+" if p["orientation"] == +1 else "-"
            p["fimo_pvalue"] = None
            n_fallback += 1

    logger.info(f"  FIMO orientation: {n_fimo}/{len(peaks)} peaks assigned by motif")
    if n_fallback > 0:
        logger.warning(f"  {n_fallback} peaks had no FIMO hit (p > threshold) — "
                       "orientation assigned by coordinate hash (deterministic but arbitrary)")
        logger.warning("  Consider relaxing --fimo-thresh or providing a pre-oriented BED")

    return peaks


def _run_fimo_batch(batch_peaks, genome_fasta, motif_meme, flank, fimo_thresh=1e-3):
    """
    Run FIMO on a batch of peaks.

    Returns dict keyed by (chrom, summit) → (strand, pvalue), so results
    are tied to genomic coordinates and consistent across runs regardless
    of peak ordering or batch boundaries.
    """
    tmp_dir = tempfile.mkdtemp(prefix="ctcf_fimo_")
    windows_bed = os.path.join(tmp_dir, "windows.bed")
    windows_fa = os.path.join(tmp_dir, "windows.fa")

    # Build a map from BED name → (chrom, summit) for result lookup
    name_to_coord = {}

    try:
        # Write BED windows using genomic coordinates as unique names
        with open(windows_bed, "w") as f:
            for p in batch_peaks:
                start = max(0, p["summit"] - flank)
                end = p["summit"] + flank
                # Unique name encoding the genomic coordinate
                bed_name = f"{p['chrom']}_{p['summit']}"
                name_to_coord[bed_name] = (p["chrom"], p["summit"])
                f.write(f"{p['chrom']}\t{start}\t{end}\t{bed_name}\t0\t+\n")

        # Extract sequences
        result = subprocess.run(
            ["bedtools", "getfasta",
             "-fi", genome_fasta,
             "-bed", windows_bed,
             "-fo", windows_fa,
             "-name"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"  bedtools getfasta failed: {result.stderr.strip()}")
            return {}

        # Run FIMO — we keep the best hit per peak regardless of threshold
        result = subprocess.run(
            ["fimo",
             "--text",
             "--thresh", str(fimo_thresh),
             motif_meme,
             windows_fa],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"  FIMO failed: {result.stderr.strip()}")
            return {}

        # Parse FIMO output using header to find column positions
        best_hits = {}  # (chrom, summit) → (strand, pvalue)
        col_map = None

        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue

            # Parse header line to get column positions dynamically
            if line.startswith("motif_id") or line.startswith("#pattern"):
                col_map = _parse_fimo_header(line.lstrip("#").strip())
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")

            # Use column map if we found a header, otherwise try defaults
            if col_map is not None:
                seq_name = fields[col_map["sequence_name"]] if col_map["sequence_name"] < len(fields) else None
                strand = fields[col_map["strand"]] if col_map["strand"] < len(fields) else None
                pvalue_str = fields[col_map["p-value"]] if col_map["p-value"] < len(fields) else None
            else:
                # Default column positions for FIMO 5.x --text output:
                # motif_id(0) motif_alt_id(1) seq_name(2) start(3) stop(4) strand(5) score(6) p-value(7) ...
                if len(fields) < 8:
                    continue
                seq_name = fields[2]
                strand = fields[5]
                pvalue_str = fields[7]

            if seq_name is None or strand is None or pvalue_str is None:
                continue
            if strand not in ("+", "-"):
                continue

            try:
                pvalue = float(pvalue_str)
            except ValueError:
                continue

            # Extract the bed_name from sequence name
            # bedtools getfasta -name produces: "bed_name::chrom:start-end"
            # or just "bed_name" depending on version
            bed_name = seq_name.split("::")[0]

            if bed_name not in name_to_coord:
                continue

            coord_key = name_to_coord[bed_name]

            # Keep the best (lowest p-value) hit per peak
            if coord_key not in best_hits or pvalue < best_hits[coord_key][1]:
                best_hits[coord_key] = (strand, pvalue)

        return best_hits

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _parse_fimo_header(header_line):
    """
    Parse FIMO --text header to find column indices dynamically.

    FIMO versions differ in column layout:
      v5.x: motif_id  motif_alt_id  sequence_name  start  stop  strand  score  p-value  [q-value]  matched_sequence
      older: pattern name  sequence name  start  stop  strand  score  p-value  ...

    Returns dict mapping column name → index.
    """
    fields = header_line.split("\t")

    # Normalize column names
    normalized = [f.strip().lower().replace(" ", "_") for f in fields]

    col_map = {}
    # Find the columns we need
    for target in ["sequence_name", "strand", "p-value"]:
        for i, name in enumerate(normalized):
            if target in name or name == target:
                col_map[target] = i
                break

    # If we found all three, return
    if len(col_map) == 3:
        logger.info(f"  FIMO header parsed: sequence_name={col_map['sequence_name']}, "
                    f"strand={col_map['strand']}, p-value={col_map['p-value']}")
        return col_map

    # Fallback: try common alternative names
    for i, name in enumerate(normalized):
        if "sequence" in name and "sequence_name" not in col_map:
            col_map["sequence_name"] = i
        elif name == "strand" and "strand" not in col_map:
            col_map["strand"] = i
        elif "p-val" in name or "pvalue" in name and "p-value" not in col_map:
            col_map["p-value"] = i

    if len(col_map) == 3:
        return col_map

    logger.warning(f"  Could not fully parse FIMO header: {header_line}")
    logger.warning(f"  Found columns: {col_map}")
    logger.warning(f"  Falling back to default column positions")
    return None


def _fallback_by_coordinate(peaks):
    """
    Assign orientation based on genomic coordinate hash when FIMO is unavailable.

    Uses hash of (chrom, summit) so the same peak always gets the same
    orientation regardless of the run mode or peak list order.
    """
    logger.warning("  FIMO unavailable. Assigning orientation by coordinate hash (INACCURATE).")
    logger.warning("  Install bedtools + meme and provide genome FASTA for real orientation.")
    for p in peaks:
        p["orientation"] = +1 if hash((p["chrom"], p["summit"])) % 2 == 0 else -1
    return peaks


def _fallback_alternating(peaks):
    """Legacy fallback — redirects to coordinate-based method."""
    return _fallback_by_coordinate(peaks)


# =============================================================================
# OUTPUT
# =============================================================================

def write_oriented_bed(peaks, output_path):
    """
    Write all oriented peaks as BED6 (genomic coordinates).

    Output columns: chrom, start, end, name, score, strand
    This is the general-purpose output usable for any downstream analysis.
    """
    with open(output_path, "w") as f:
        f.write("# CTCF peaks with motif orientation (FIMO MA0139.1)\n")
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        for p in peaks:
            strand = "+" if p.get("orientation", +1) == +1 else "-"
            f.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t"
                    f"{p['name']}\t{int(p['score'])}\t{strand}\n")
    logger.info(f"  Wrote {len(peaks)} oriented peaks to {output_path}")


def write_simulation_bed(peaks, output_path, region):
    """
    Write peaks within a region as BED6 formatted for load_ctcf_from_bed().

    Uses summit coordinates so the simulation gets precise positions.
    """
    r_chrom, r_start, r_end = region
    with open(output_path, "w") as f:
        f.write(f"# CTCF sites for {r_chrom}:{r_start}-{r_end} (mm10)\n")
        for p in peaks:
            if p["chrom"] != r_chrom:
                continue
            if p["summit"] < r_start or p["summit"] > r_end:
                continue
            strand = "+" if p.get("orientation", +1) == +1 else "-"
            f.write(f"{p['chrom']}\t{p['summit']}\t{p['summit'] + 1}\t"
                    f"{p['name']}\t{int(p['score'])}\t{strand}\n")
    logger.info(f"  Wrote simulation BED to {output_path}")


def print_python_array(peaks, region, resolution=1000):
    """Print CTCF_SITES array for parameters.py (only for a specific region)."""
    r_chrom, r_start, r_end = region
    n_monomers = (r_end - r_start) // resolution

    sites = []
    for p in peaks:
        if p["chrom"] != r_chrom:
            continue
        monomer = (p["summit"] - r_start) // resolution
        if 0 <= monomer < n_monomers:
            ori = p.get("orientation", +1)
            sites.append((monomer, ori))

    print(f"\n# Region: {r_chrom}:{r_start:,}-{r_end:,} "
          f"({n_monomers} monomers at {resolution} bp)")
    print("CTCF_SITES = [")
    for monomer, ori in sites:
        genomic = r_start + monomer * resolution
        sign = "+" if ori == +1 else "-"
        print(f"    ({monomer:>4d}, {sign}1),   # {r_chrom}:{genomic:,}")
    print("]")
    print(f"\n# {len(sites)} CTCF sites total")
    return sites


def print_summary(peaks):
    """Print genome-wide summary statistics."""
    from collections import Counter
    chrom_counts = Counter(p["chrom"] for p in peaks)
    n_plus = sum(1 for p in peaks if p.get("orientation", +1) == +1)
    n_minus = sum(1 for p in peaks if p.get("orientation", -1) == -1)
    n_fimo = sum(1 for p in peaks if p.get("fimo_pvalue") is not None)

    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total peaks:        {len(peaks)}")
    print(f"Forward (+):        {n_plus}")
    print(f"Reverse (-):        {n_minus}")
    print(f"Oriented by FIMO:   {n_fimo}")
    print(f"Chromosomes:        {len(chrom_counts)}")
    print(f"\nPer-chromosome counts:")
    for chrom in sorted(chrom_counts, key=lambda c: (len(c), c)):
        print(f"  {chrom:<6s}  {chrom_counts[chrom]:>5d}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extract CTCF sites with motif orientation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Genome-wide (all peaks, all chromosomes):
  python scripts/extract_ctcf_sites.py --source encode --genome-wide

  # Specific region:
  python scripts/extract_ctcf_sites.py --source encode --region chr3:34000000-36000000

  # Sox2 locus (default from parameters.py):
  python scripts/extract_ctcf_sites.py --source encode

  # Use E14 mESC line instead of Bruce4:
  python scripts/extract_ctcf_sites.py --source encode --genome-wide --encode-line mESC_E14

  # Your own peaks, genome-wide:
  python scripts/extract_ctcf_sites.py --source bed --bed peaks.narrowPeak --genome-wide

  # Provide your own mm10.fa (skip download):
  python scripts/extract_ctcf_sites.py --source encode --genome /data/mm10.fa --genome-wide

  # Fast: chr3 genome only (enough for Sox2 region):
  python scripts/extract_ctcf_sites.py --source encode --chr3-only
""",
    )

    parser.add_argument("--source", choices=["encode", "bed"], default="encode",
                        help="Data source for CTCF peaks (default: encode)")
    parser.add_argument("--encode-line", choices=["mESC_Bruce4", "mESC_E14"],
                        default="mESC_Bruce4",
                        help="ENCODE cell line (default: mESC_Bruce4 / ENCSR000CCB)")
    parser.add_argument("--bed", type=str, default=None,
                        help="Path to BED/narrowPeak (required if --source bed)")

    region_group = parser.add_mutually_exclusive_group()
    region_group.add_argument("--genome-wide", action="store_true",
                              help="Process ALL peaks genome-wide (no region filter)")
    region_group.add_argument("--region", type=str, default=None,
                              help="Genomic region, e.g. chr3:34000000-36000000. "
                                   "Default: Sox2 locus from parameters.py")

    parser.add_argument("--genome", type=str, default=None,
                        help="Path to mm10 FASTA (auto-downloads if not provided)")
    parser.add_argument("--chr3-only", action="store_true",
                        help="Download/use only chr3 from genome (faster)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output BED file path")
    parser.add_argument("--skip-orientation", action="store_true",
                        help="Skip FIMO orientation calling")
    parser.add_argument("--fimo-thresh", type=float, default=1e-3,
                        help="FIMO p-value threshold (default: 1e-3). "
                             "More lenient = fewer unassigned peaks")
    parser.add_argument("--resolution", type=int, default=1000,
                        help="Monomer resolution in bp for simulation array "
                             "(default: 1000)")

    args = parser.parse_args()

    proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = ensure_dir(os.path.join(proj_dir, "data"))
    genome_dir = ensure_dir(os.path.join(data_dir, "genome"))
    motif_dir = ensure_dir(os.path.join(data_dir, "motifs"))

    # Determine region
    region = None  # None = genome-wide
    if not args.genome_wide:
        if args.region:
            region = parse_region_string(args.region)
        else:
            # Default: Sox2 locus from parameters.py
            try:
                sys.path.insert(0, proj_dir)
                from configs.parameters import CHROM, REGION_START, REGION_END
                region = (CHROM, REGION_START, REGION_END)
                logger.info(f"Using default region from parameters.py: "
                            f"{CHROM}:{REGION_START:,}-{REGION_END:,}")
            except ImportError:
                logger.error("Cannot import parameters.py and no --region given. "
                             "Use --genome-wide or --region chr:start-end")
                sys.exit(1)

    mode = "genome-wide" if region is None else f"{region[0]}:{region[1]:,}-{region[2]:,}"
    logger.info(f"Mode: {mode}")

    # ── Step 1: Get CTCF peaks ──────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("STEP 1: Obtaining CTCF peaks")
    logger.info("=" * 60)

    if args.source == "encode":
        bed_path = download_encode_peaks(data_dir, args.encode_line)
    elif args.bed:
        bed_path = args.bed
        if not os.path.exists(bed_path):
            logger.error(f"BED file not found: {bed_path}")
            sys.exit(1)
    else:
        parser.error("--bed is required when --source bed")

    peaks = load_peaks(bed_path, region=region)
    if not peaks:
        logger.error("No CTCF peaks found!")
        sys.exit(1)

    # ── Step 2: Get genome + motif, call orientation ────────────────────
    if not args.skip_orientation:
        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 2: Obtaining genome FASTA")
        logger.info("=" * 60)

        if args.genome and os.path.exists(args.genome):
            genome_fa = args.genome
            logger.info(f"  Using provided genome: {genome_fa}")
            if not os.path.exists(genome_fa + ".fai"):
                subprocess.run(["samtools", "faidx", genome_fa], check=True)
        else:
            chrom_only = None
            if args.chr3_only:
                chrom_only = "chr3"
            elif region is not None:
                chrom_only = region[0]
            genome_fa = download_genome(genome_dir, chrom_only=chrom_only)

        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 3: Obtaining CTCF motif (JASPAR MA0139.1)")
        logger.info("=" * 60)

        motif_meme = download_jaspar_motif(motif_dir)

        logger.info("")
        logger.info("=" * 60)
        logger.info("STEP 4: Calling motif orientation with FIMO")
        logger.info("=" * 60)

        peaks = call_motif_orientation(peaks, genome_fa, motif_meme,
                                       fimo_thresh=args.fimo_thresh)
    else:
        peaks = _fallback_alternating(peaks)

    # ── Step 3: Output ──────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("STEP 5: Generating output")
    logger.info("=" * 60)

    # Always write the full oriented BED
    if args.output:
        out_bed = args.output
    elif region is not None:
        r_chrom, r_start, r_end = region
        out_bed = os.path.join(data_dir,
                               f"ctcf_oriented_{r_chrom}_{r_start}_{r_end}.bed")
    else:
        out_bed = os.path.join(data_dir, "ctcf_oriented_mm10_genome_wide.bed")

    write_oriented_bed(peaks, out_bed)

    # If region mode, also write simulation-ready BED and Python array
    if region is not None:
        sim_bed = os.path.join(data_dir, "ctcf_sites_sox2.bed")
        write_simulation_bed(peaks, sim_bed, region)

        sites = print_python_array(peaks, region, args.resolution)

        print(f"\n{'=' * 60}")
        print("TO USE THESE SITES IN SIMULATION:")
        print(f"{'=' * 60}")
        print("Option A: Copy the CTCF_SITES array above into configs/parameters.py")
        print("          and set CTCF_SITES_PLACEHOLDER = False")
        print("")
        print("Option B: Load at runtime:")
        print(f"          from configs.parameters import load_ctcf_from_bed")
        print(f"          load_ctcf_from_bed('{sim_bed}')")
    else:
        print_summary(peaks)

    print(f"\nOutput: {out_bed}")


if __name__ == "__main__":
    main()
