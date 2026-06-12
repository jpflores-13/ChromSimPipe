#!/usr/bin/env python
"""
Extract CTCF binding sites with motif orientation from ENCODE ChIP-seq data.
Human genome (GRCh38/hg38).

Can run genome-wide or on a specific region. When a region is specified,
also outputs a CTCF_SITES array for polychrom simulations.

Data sources:
  - ENCODE CTCF ChIP-seq (GRCh38), multiple human cell lines:
    GM12878:  ENCSR000AKB → ENCFF356LIU (optimal IDR peaks)
    K562:     ENCSR000AKO → ENCFF119XFJ
    H1-hESC:              → ENCFF368LWM
    HCT116:   ENCSR240PRQ  (Rao et al. 2017 line)
    IMR-90:   ENCSR000EFI
    HeLa-S3:  ENCSR000AOA

  - CTCF motif orientation via FIMO:
    JASPAR motif MA0139.1 (auto-downloaded)
    hg38 genome FASTA (auto-downloaded from UCSC)

Environment setup:
    mamba env create -f envs/ctcf_extraction.yml
    conda activate ctcf_extraction

Usage:
    # Genome-wide, GM12878 (default):
    python scripts/extract_ctcf_sites_hg38.py --source encode --genome-wide

    # Specific region:
    python scripts/extract_ctcf_sites_hg38.py --source encode --region chr1:50000000-52000000

    # Different cell line:
    python scripts/extract_ctcf_sites_hg38.py --source encode --cell-line K562 --genome-wide

    # Your own peaks:
    python scripts/extract_ctcf_sites_hg38.py --source bed --bed my_peaks.narrowPeak --genome-wide
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
# CONSTANTS & URLS — hg38
# =============================================================================

GENOME_ASSEMBLY = "hg38"

ENCODE_CTCF_PEAKS = {
    "GM12878": {
        "experiment": "ENCSR000AKB",
        "file": "ENCFF356LIU",
        "url": "https://www.encodeproject.org/files/ENCFF356LIU/@@download/ENCFF356LIU.bed.gz",
        "description": "CTCF ChIP-seq optimal IDR peaks, GM12878, GRCh38",
    },
    "K562": {
        "experiment": "ENCSR000AKO",
        "file": "ENCFF119XFJ",
        "url": "https://www.encodeproject.org/files/ENCFF119XFJ/@@download/ENCFF119XFJ.bed.gz",
        "description": "CTCF ChIP-seq optimal IDR peaks, K562, GRCh38",
    },
    "H1-hESC": {
        "experiment": "ENCSR000ATN",
        "file": "ENCFF368LWM",
        "url": "https://www.encodeproject.org/files/ENCFF368LWM/@@download/ENCFF368LWM.bed.gz",
        "description": "CTCF ChIP-seq optimal IDR peaks, H1-hESC, GRCh38",
    },
    "HCT116": {
        "experiment": "ENCSR240PRQ",
        "file": "ENCSR240PRQ",  # exact file TBD — download experiment page
        "url": None,  # see note below
        "description": "CTCF ChIP-seq, HCT116, GRCh38 (Broad/Bernstein)",
    },
    "IMR-90": {
        "experiment": "ENCSR000EFI",
        "file": "ENCSR000EFI",
        "url": None,
        "description": "CTCF ChIP-seq, IMR-90, GRCh38 (Stanford/Snyder)",
    },
    "HeLa-S3": {
        "experiment": "ENCSR000AOA",
        "file": "ENCSR000AOA",
        "url": None,
        "description": "CTCF ChIP-seq, HeLa-S3, GRCh38 (Broad)",
    },
}

# Lines with confirmed direct-download narrowPeak URLs on GRCh38
AVAILABLE_CELL_LINES = ["GM12878", "K562", "H1-hESC"]

JASPAR_CTCF_MEME_URL = (
    "https://jaspar.elixir.no/api/v1/matrix/MA0139.1/?format=meme"
)

HG38_FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"


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


def download_encode_peaks(output_dir, cell_type="GM12878"):
    """Download CTCF narrowPeak from ENCODE and decompress."""
    if cell_type not in ENCODE_CTCF_PEAKS:
        raise ValueError(f"Unknown cell type '{cell_type}'. "
                         f"Available: {list(ENCODE_CTCF_PEAKS.keys())}")

    info = ENCODE_CTCF_PEAKS[cell_type]

    if info["url"] is None:
        logger.error(f"  No direct download URL for {cell_type}.")
        logger.error(f"  Visit https://www.encodeproject.org/experiments/{info['experiment']}/")
        logger.error(f"  Download the optimal IDR narrowPeak (GRCh38) manually,")
        logger.error(f"  then use --source bed --bed /path/to/file.bed.gz")
        sys.exit(1)

    out_gz = os.path.join(output_dir, f"CTCF_{cell_type}_{GENOME_ASSEMBLY}_peaks.bed.gz")
    out_bed = os.path.join(output_dir, f"CTCF_{cell_type}_{GENOME_ASSEMBLY}_peaks.bed")

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
    """Download hg38 genome FASTA and index it."""
    genome_dir = ensure_dir(output_dir)
    fasta_path = os.path.join(genome_dir, "hg38.fa")
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

        download_file(HG38_FASTA_URL, gz_path,
                      "hg38 genome FASTA (~900 MB compressed)")

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
    """Parse "chr1:50000000-52000000" into ("chr1", 50000000, 52000000)."""
    region_str = region_str.replace(",", "")
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def load_peaks(bed_path, region=None):
    """Load peaks from a narrowPeak/BED, optionally filtered to a region."""
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

            if region is not None:
                r_chrom, r_start, r_end = region
                if chrom != r_chrom:
                    continue
                if end < r_start or start > r_end:
                    continue

            name = fields[3] if len(fields) > 3 else "."
            score = float(fields[4]) if len(fields) > 4 else 0
            strand = fields[5] if len(fields) > 5 else "."

            if len(fields) >= 10:
                try:
                    summit_offset = int(fields[9])
                    summit = start + summit_offset
                except ValueError:
                    summit = (start + end) // 2
            else:
                summit = (start + end) // 2

            peaks.append({
                "chrom": chrom, "start": start, "end": end, "summit": summit,
                "name": name, "score": score, "strand": strand,
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
# MOTIF ORIENTATION (shared logic with mm10 script)
# =============================================================================

def call_motif_orientation(peaks, genome_fasta, motif_meme, flank=100, fimo_thresh=1e-3):
    """Determine CTCF motif orientation at each peak using FIMO."""
    if not peaks:
        return peaks

    has_strand = all(p["strand"] in ("+", "-") for p in peaks)
    if has_strand:
        logger.info("  BED file already has strand info — using directly")
        for p in peaks:
            p["orientation"] = +1 if p["strand"] == "+" else -1
        return peaks

    for tool in ["bedtools", "fimo"]:
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            logger.error(f"  {tool} not found!")
            return _fallback_by_coordinate(peaks)

    logger.info(f"  Calling motif orientation with FIMO on {len(peaks)} peaks "
                f"(±{flank} bp windows)...")

    BATCH_SIZE = 50_000
    all_best_hits = {}

    for batch_start in range(0, len(peaks), BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, len(peaks))
        batch = peaks[batch_start:batch_end]

        if len(peaks) > BATCH_SIZE:
            logger.info(f"  Processing batch {batch_start // BATCH_SIZE + 1} "
                        f"({batch_start}-{batch_end} of {len(peaks)})...")

        batch_hits = _run_fimo_batch(batch, genome_fasta, motif_meme, flank,
                                     fimo_thresh)
        all_best_hits.update(batch_hits)

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
            p["orientation"] = +1 if hash((p["chrom"], p["summit"])) % 2 == 0 else -1
            p["strand"] = "+" if p["orientation"] == +1 else "-"
            p["fimo_pvalue"] = None
            n_fallback += 1

    logger.info(f"  FIMO orientation: {n_fimo}/{len(peaks)} peaks assigned by motif")
    if n_fallback > 0:
        logger.warning(f"  {n_fallback} peaks had no FIMO hit — "
                       "orientation assigned by coordinate hash")

    return peaks


def _run_fimo_batch(batch_peaks, genome_fasta, motif_meme, flank, fimo_thresh=1e-3):
    """Run FIMO on a batch. Returns {(chrom, summit): (strand, pvalue)}."""
    tmp_dir = tempfile.mkdtemp(prefix="ctcf_fimo_hg38_")
    windows_bed = os.path.join(tmp_dir, "windows.bed")
    windows_fa = os.path.join(tmp_dir, "windows.fa")

    name_to_coord = {}

    try:
        with open(windows_bed, "w") as f:
            for p in batch_peaks:
                start = max(0, p["summit"] - flank)
                end = p["summit"] + flank
                bed_name = f"{p['chrom']}_{p['summit']}"
                name_to_coord[bed_name] = (p["chrom"], p["summit"])
                f.write(f"{p['chrom']}\t{start}\t{end}\t{bed_name}\t0\t+\n")

        result = subprocess.run(
            ["bedtools", "getfasta", "-fi", genome_fasta,
             "-bed", windows_bed, "-fo", windows_fa, "-name"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"  bedtools getfasta failed: {result.stderr.strip()}")
            return {}

        result = subprocess.run(
            ["fimo", "--text", "--thresh", str(fimo_thresh),
             motif_meme, windows_fa],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            logger.error(f"  FIMO failed: {result.stderr.strip()}")
            return {}

        best_hits = {}
        col_map = None

        for line in result.stdout.strip().split("\n"):
            if not line.strip():
                continue
            if line.startswith("motif_id") or line.startswith("#pattern"):
                col_map = _parse_fimo_header(line.lstrip("#").strip())
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if col_map is not None:
                seq_name = fields[col_map["sequence_name"]] if col_map["sequence_name"] < len(fields) else None
                strand = fields[col_map["strand"]] if col_map["strand"] < len(fields) else None
                pvalue_str = fields[col_map["p-value"]] if col_map["p-value"] < len(fields) else None
            else:
                if len(fields) < 8:
                    continue
                seq_name, strand, pvalue_str = fields[2], fields[5], fields[7]

            if not all([seq_name, strand, pvalue_str]) or strand not in ("+", "-"):
                continue

            try:
                pvalue = float(pvalue_str)
            except ValueError:
                continue

            bed_name = seq_name.split("::")[0]
            if bed_name not in name_to_coord:
                continue

            coord_key = name_to_coord[bed_name]
            if coord_key not in best_hits or pvalue < best_hits[coord_key][1]:
                best_hits[coord_key] = (strand, pvalue)

        return best_hits

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _parse_fimo_header(header_line):
    """Parse FIMO header to find column indices dynamically."""
    fields = header_line.split("\t")
    normalized = [f.strip().lower().replace(" ", "_") for f in fields]

    col_map = {}
    for target in ["sequence_name", "strand", "p-value"]:
        for i, name in enumerate(normalized):
            if target in name or name == target:
                col_map[target] = i
                break

    if len(col_map) == 3:
        logger.info(f"  FIMO header parsed: {col_map}")
        return col_map

    for i, name in enumerate(normalized):
        if "sequence" in name and "sequence_name" not in col_map:
            col_map["sequence_name"] = i
        elif name == "strand" and "strand" not in col_map:
            col_map["strand"] = i
        elif ("p-val" in name or "pvalue" in name) and "p-value" not in col_map:
            col_map["p-value"] = i

    if len(col_map) == 3:
        return col_map

    logger.warning(f"  Could not parse FIMO header: {header_line}")
    return None


def _fallback_by_coordinate(peaks):
    """Deterministic fallback based on genomic coordinate hash."""
    logger.warning("  FIMO unavailable. Assigning orientation by coordinate hash (INACCURATE).")
    for p in peaks:
        p["orientation"] = +1 if hash((p["chrom"], p["summit"])) % 2 == 0 else -1
    return peaks


# =============================================================================
# OUTPUT
# =============================================================================

def write_oriented_bed(peaks, output_path):
    """Write all oriented peaks as BED6."""
    with open(output_path, "w") as f:
        f.write(f"# CTCF peaks with motif orientation (FIMO MA0139.1) — {GENOME_ASSEMBLY}\n")
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        for p in peaks:
            strand = "+" if p.get("orientation", +1) == +1 else "-"
            f.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t"
                    f"{p['name']}\t{int(p['score'])}\t{strand}\n")
    logger.info(f"  Wrote {len(peaks)} oriented peaks to {output_path}")


def write_simulation_bed(peaks, output_path, region):
    """Write peaks within a region as BED6 for load_ctcf_from_bed()."""
    r_chrom, r_start, r_end = region
    with open(output_path, "w") as f:
        f.write(f"# CTCF sites for {r_chrom}:{r_start}-{r_end} ({GENOME_ASSEMBLY})\n")
        for p in peaks:
            if p["chrom"] != r_chrom or p["summit"] < r_start or p["summit"] > r_end:
                continue
            strand = "+" if p.get("orientation", +1) == +1 else "-"
            f.write(f"{p['chrom']}\t{p['summit']}\t{p['summit'] + 1}\t"
                    f"{p['name']}\t{int(p['score'])}\t{strand}\n")
    logger.info(f"  Wrote simulation BED to {output_path}")


def print_python_array(peaks, region, resolution=1000):
    """Print CTCF_SITES array for a specific region."""
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

    print(f"\n# Region: {r_chrom}:{r_start:,}-{r_end:,} ({GENOME_ASSEMBLY})")
    print(f"# {n_monomers} monomers at {resolution} bp")
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
    print(f"SUMMARY ({GENOME_ASSEMBLY})")
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
        description=f"Extract CTCF sites with motif orientation ({GENOME_ASSEMBLY})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  # Genome-wide, GM12878 (default):
  python scripts/extract_ctcf_sites_hg38.py --source encode --genome-wide

  # Specific region:
  python scripts/extract_ctcf_sites_hg38.py --source encode --region chr1:50000000-52000000

  # K562 cell line:
  python scripts/extract_ctcf_sites_hg38.py --source encode --cell-line K562 --genome-wide

  # H1-hESC:
  python scripts/extract_ctcf_sites_hg38.py --source encode --cell-line H1-hESC --genome-wide

  # Your own peaks:
  python scripts/extract_ctcf_sites_hg38.py --source bed --bed peaks.narrowPeak --genome-wide

  # Provide hg38.fa (skip download):
  python scripts/extract_ctcf_sites_hg38.py --source encode --genome /data/hg38.fa --genome-wide

Available ENCODE cell lines with direct download:
  GM12878  (ENCSR000AKB / ENCFF356LIU) — lymphoblastoid
  K562     (ENCSR000AKO / ENCFF119XFJ) — CML
  H1-hESC  (ENCSR000ATN / ENCFF368LWM) — embryonic stem cell

Other experiments (download narrowPeak manually, use --source bed):
  HCT116   ENCSR240PRQ
  IMR-90   ENCSR000EFI
  HeLa-S3  ENCSR000AOA
""",
    )

    parser.add_argument("--source", choices=["encode", "bed"], default="encode")
    parser.add_argument("--cell-line",
                        choices=list(ENCODE_CTCF_PEAKS.keys()),
                        default="GM12878",
                        help=f"ENCODE cell line (default: GM12878). "
                             f"Direct download: {', '.join(AVAILABLE_CELL_LINES)}")
    parser.add_argument("--bed", type=str, default=None,
                        help="Path to BED/narrowPeak (required if --source bed)")

    region_group = parser.add_mutually_exclusive_group()
    region_group.add_argument("--genome-wide", action="store_true",
                              help="Process ALL peaks genome-wide")
    region_group.add_argument("--region", type=str, default=None,
                              help="Genomic region, e.g. chr1:50000000-52000000")

    parser.add_argument("--genome", type=str, default=None,
                        help="Path to hg38 FASTA (auto-downloads if not provided)")
    parser.add_argument("--chrom-only", type=str, default=None,
                        help="Download/use only this chromosome (e.g. chr1)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output BED file path")
    parser.add_argument("--skip-orientation", action="store_true",
                        help="Skip FIMO orientation calling")
    parser.add_argument("--fimo-thresh", type=float, default=1e-3,
                        help="FIMO p-value threshold (default: 1e-3)")
    parser.add_argument("--resolution", type=int, default=1000,
                        help="Monomer resolution in bp (default: 1000)")

    args = parser.parse_args()

    proj_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = ensure_dir(os.path.join(proj_dir, "data"))
    genome_dir = ensure_dir(os.path.join(data_dir, "genome"))
    motif_dir = ensure_dir(os.path.join(data_dir, "motifs"))

    # Determine region
    region = None
    if not args.genome_wide:
        if args.region:
            region = parse_region_string(args.region)
        else:
            logger.error("Specify --region chr:start-end or --genome-wide")
            sys.exit(1)

    mode = "genome-wide" if region is None else f"{region[0]}:{region[1]:,}-{region[2]:,}"
    logger.info(f"Assembly: {GENOME_ASSEMBLY}")
    logger.info(f"Mode: {mode}")

    # ── Step 1: Get CTCF peaks ──────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("STEP 1: Obtaining CTCF peaks")
    logger.info("=" * 60)

    if args.source == "encode":
        bed_path = download_encode_peaks(data_dir, args.cell_line)
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

    # ── Step 2: Genome + FIMO ───────────────────────────────────────────
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
            chrom_only = args.chrom_only
            if chrom_only is None and region is not None:
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
        peaks = _fallback_by_coordinate(peaks)

    # ── Step 3: Output ──────────────────────────────────────────────────
    logger.info("")
    logger.info("=" * 60)
    logger.info("STEP 5: Generating output")
    logger.info("=" * 60)

    if args.output:
        out_bed = args.output
    elif region is not None:
        r_chrom, r_start, r_end = region
        out_bed = os.path.join(data_dir,
                               f"ctcf_oriented_{GENOME_ASSEMBLY}_{r_chrom}_{r_start}_{r_end}.bed")
    else:
        cell = args.cell_line if args.source == "encode" else "custom"
        out_bed = os.path.join(data_dir,
                               f"ctcf_oriented_{GENOME_ASSEMBLY}_{cell}_genome_wide.bed")

    write_oriented_bed(peaks, out_bed)

    if region is not None:
        r_chrom, r_start, r_end = region
        sim_bed = os.path.join(data_dir,
                               f"ctcf_sites_{GENOME_ASSEMBLY}_{r_chrom}_{r_start}_{r_end}.bed")
        write_simulation_bed(peaks, sim_bed, region)
        sites = print_python_array(peaks, region, args.resolution)

        print(f"\n{'=' * 60}")
        print("TO USE THESE SITES IN SIMULATION:")
        print(f"{'=' * 60}")
        print("Option A: Copy the CTCF_SITES array above into your parameters file")
        print("")
        print("Option B: Load at runtime:")
        print(f"          load_ctcf_from_bed('{sim_bed}')")
    else:
        print_summary(peaks)

    print(f"\nOutput: {out_bed}")


if __name__ == "__main__":
    main()
