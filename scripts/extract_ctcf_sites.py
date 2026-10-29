#!/usr/bin/env python
"""
Orient CTCF peaks by motif strand for loop-extrusion simulations.

Reads CTCF peaks (narrowPeak/BED) inside a region, scans ± --flank bp around
each summit with FIMO for the CTCF motif in --motif (MEME format), and writes
a BED6 whose strand column is the strand of the best motif hit. Peaks that
already carry a +/- strand are used as given.

Called by the orientCTCF rule in workflows/ChromSimPipe.snakefile; every
argument comes from config/ChromSimConfig.yaml and the samplesheet.

Usage:
    python scripts/extract_ctcf_sites.py \
        --bed peaks.narrowPeak --region chr1:64210000-66170000 \
        --genome genome.fa --motif MA0139.1.meme --assembly hg38 \
        --output ctcf_oriented.bed
"""

import os
import sys
import argparse
import logging
import shutil
import subprocess
import tempfile

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


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
# MOTIF ORIENTATION
# =============================================================================

def call_motif_orientation(peaks, genome_fasta, motif_meme, flank=100,
                           fimo_thresh=1e-3, unoriented="drop"):
    """
    Determine CTCF motif orientation at each peak using FIMO.

    Peaks that already carry a +/- strand are used as given. Peaks with no
    FIMO hit have no measured orientation: they are dropped
    (``unoriented='drop'``) or the run stops (``unoriented='error'``).
    """
    if not peaks:
        return peaks

    has_strand = all(p["strand"] in ("+", "-") for p in peaks)
    if has_strand:
        logger.info("  BED file already has strand info — using directly")
        for p in peaks:
            p["orientation"] = +1 if p["strand"] == "+" else -1
        return peaks

    for tool in ["bedtools", "fimo"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} not found on PATH — check the ctcf conda env")

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

        all_best_hits.update(
            _run_fimo_batch(batch, genome_fasta, motif_meme, flank, fimo_thresh)
        )

    oriented, missing = [], []
    for p in peaks:
        key = (p["chrom"], p["summit"])
        if key in all_best_hits:
            strand, pval = all_best_hits[key]
            p["orientation"] = +1 if strand == "+" else -1
            p["strand"] = strand
            p["fimo_pvalue"] = pval
            oriented.append(p)
        else:
            missing.append(p)

    logger.info(f"  FIMO orientation: {len(oriented)}/{len(peaks)} peaks assigned by motif")
    if missing:
        listing = ", ".join(f"{p['chrom']}:{p['summit']}" for p in missing)
        if unoriented == "error":
            raise RuntimeError(
                f"{len(missing)} peaks had no FIMO hit (p < {fimo_thresh}): {listing}. "
                "Set ctcf_motif.unoriented: 'drop' to leave them out."
            )
        logger.warning(f"  Dropped {len(missing)} peaks with no FIMO hit "
                       f"(p < {fimo_thresh}): {listing}")

    return oriented


def _run_fimo_batch(batch_peaks, genome_fasta, motif_meme, flank, fimo_thresh=1e-3):
    """Run FIMO on a batch. Returns {(chrom, summit): (strand, pvalue)}."""
    tmp_dir = tempfile.mkdtemp(prefix="ctcf_fimo_")
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
            raise RuntimeError(f"bedtools getfasta failed: {result.stderr.strip()}")

        result = subprocess.run(
            ["fimo", "--text", "--thresh", str(fimo_thresh),
             motif_meme, windows_fa],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"FIMO failed: {result.stderr.strip()}")

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


# =============================================================================
# OUTPUT
# =============================================================================

def write_oriented_bed(peaks, output_path, assembly, motif_name):
    """Write all oriented peaks as BED6."""
    with open(output_path, "w") as f:
        f.write(f"# CTCF peaks with motif orientation (FIMO {motif_name}) — {assembly}\n")
        f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
        for p in peaks:
            strand = "+" if p.get("orientation", +1) == +1 else "-"
            f.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t"
                    f"{p['name']}\t{int(p['score'])}\t{strand}\n")
    logger.info(f"  Wrote {len(peaks)} oriented peaks to {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def _motif_name(meme_path):
    """First MOTIF id in a MEME file (for the output header)."""
    with open(meme_path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                return " ".join(line.split()[1:])
    raise ValueError(f"No MOTIF entry in {meme_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Orient CTCF peaks by FIMO motif strand",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bed", required=True, help="CTCF peaks (narrowPeak/BED)")
    parser.add_argument("--region", required=True, help="chrom:start-end")
    parser.add_argument("--genome", required=True, help="Genome FASTA")
    parser.add_argument("--motif", required=True, help="CTCF motif, MEME format")
    parser.add_argument("--assembly", required=True, help="Assembly name for the BED header")
    parser.add_argument("--output", required=True, help="Output BED6 path")
    parser.add_argument("--fimo-thresh", type=float, required=True,
                        help="FIMO p-value threshold")
    parser.add_argument("--flank", type=int, required=True,
                        help="± bp around each summit scanned by FIMO")
    parser.add_argument("--unoriented", choices=["drop", "error"], required=True,
                        help="What to do with peaks that have no FIMO hit")
    parser.add_argument("--resolution", type=int, required=True,
                        help="bp per monomer (for the site summary)")
    args = parser.parse_args()

    region = parse_region_string(args.region)
    logger.info(f"Assembly: {args.assembly}")
    logger.info(f"Region:   {region[0]}:{region[1]:,}-{region[2]:,}")

    peaks = load_peaks(args.bed, region=region)
    if not peaks:
        logger.error(f"No CTCF peaks in {args.region} in {args.bed}")
        sys.exit(1)

    if not os.path.exists(args.genome + ".fai"):
        logger.info(f"  Indexing {args.genome} with samtools faidx")
        subprocess.run(["samtools", "faidx", args.genome], check=True)

    peaks = call_motif_orientation(
        peaks, args.genome, args.motif, flank=args.flank,
        fimo_thresh=args.fimo_thresh, unoriented=args.unoriented,
    )
    if not peaks:
        logger.error("No peaks could be oriented")
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    write_oriented_bed(peaks, args.output, args.assembly, _motif_name(args.motif))

    r_chrom, r_start, r_end = region
    n_monomers = (r_end - r_start) // args.resolution
    n_plus = sum(1 for p in peaks if p["orientation"] == +1)
    logger.info(f"  {len(peaks)} oriented sites ({n_plus} +, {len(peaks) - n_plus} -) "
                f"over {n_monomers} monomers at {args.resolution} bp")


if __name__ == "__main__":
    main()
