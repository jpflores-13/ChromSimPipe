#!/usr/bin/env python
"""
Download and extract Hi-C contact maps from Bonev et al. 2017 (GSE96107).

Downloads .mcool files (or fetches from 4DN Data Portal) and extracts the
Sox2 locus region at 1 kb resolution for comparison with simulations.

Requirements:
    pip install cooler cooltools requests
"""

import os
import sys
import logging
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from configs.parameters import CHROM, REGION_START, REGION_END, RESOLUTION, DATA

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Region string for cooler
REGION = f"{CHROM}:{REGION_START}-{REGION_END}"


def download_from_4dn(output_dir: str):
    """
    Attempt to download mcool files from the 4DN Data Portal.

    The Bonev 2017 data (GSE96107) may also be available on the 4DN portal.
    If direct mcool URLs are known, they can be specified here.
    """
    import requests

    # 4DN Data Portal search for Bonev 2017 datasets
    # These are placeholder URLs — the actual mcool files need to be located
    # through the 4DN portal (https://data.4dnucleome.org/)
    logger.info("Searching 4DN Data Portal for Bonev 2017 mcool files...")
    logger.info("Visit https://data.4dnucleome.org/ and search for GSE96107 or Bonev 2017")
    logger.info("Or use the GEO accession GSE96107 to find processed mcool files")

    return None


def download_from_geo(output_dir: str):
    """
    Download processed Hi-C data from GEO (GSE96107).

    The raw data requires significant processing (alignment, filtering, binning).
    Instead, we recommend using pre-processed mcool files if available,
    or processing with the distiller pipeline.
    """
    logger.info("="*70)
    logger.info("DATA DOWNLOAD INSTRUCTIONS")
    logger.info("="*70)
    logger.info("")
    logger.info("The Bonev et al. 2017 Hi-C data (GSE96107) needs to be obtained from:")
    logger.info("")
    logger.info("Option 1: Pre-processed mcool files")
    logger.info("  - Check the Bonev lab resources: https://www.bonevlab.com/resources")
    logger.info("  - Check 4DN Data Portal: https://data.4dnucleome.org/")
    logger.info("  - Check the HiC reanalysis repo: https://github.com/lldelisle/Hi-C_reanalysis_Bonev_2017")
    logger.info("")
    logger.info("Option 2: Process from raw FASTQ")
    logger.info(f"  - Download from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={DATA['geo_accession']}")
    logger.info("  - Process with distiller-nf or HiC-Pro pipeline")
    logger.info("  - Generate mcool files with cooler")
    logger.info("")
    logger.info("Option 3: Use cooltools to fetch from a public cooler server")
    logger.info("  - Check if data is hosted on higlass.io or similar")
    logger.info("")
    logger.info("Once you have mcool files, place them in:")
    logger.info(f"  {output_dir}/")
    logger.info("  With names: ES.mcool, CN.mcool")
    logger.info("")

    return None


def extract_region_from_mcool(mcool_path: str, resolution: int = 1000) -> np.ndarray:
    """
    Extract the Sox2 locus contact matrix from a mcool file.

    Parameters
    ----------
    mcool_path : str
        Path to .mcool file.
    resolution : int
        Bin resolution in bp.

    Returns
    -------
    matrix : np.ndarray
        Contact matrix for the region.
    """
    import cooler

    uri = f"{mcool_path}::resolutions/{resolution}"
    clr = cooler.Cooler(uri)

    logger.info(f"Extracting region {REGION} at {resolution} bp resolution")
    matrix = clr.matrix(balance=True).fetch(REGION)

    # Replace NaN with 0
    matrix = np.nan_to_num(matrix, nan=0.0)

    logger.info(f"Extracted matrix shape: {matrix.shape}")
    return matrix


def extract_observed_expected(mcool_path: str, resolution: int = 1000):
    """
    Compute observed/expected matrix and P(s) curve for the region.

    Parameters
    ----------
    mcool_path : str
        Path to .mcool file.
    resolution : int
        Bin resolution.

    Returns
    -------
    obs_exp : np.ndarray
        Observed/expected contact matrix.
    ps_curve : tuple of (distances, contact_probs)
    """
    import cooler
    import cooltools

    uri = f"{mcool_path}::resolutions/{resolution}"
    clr = cooler.Cooler(uri)

    # Compute expected (distance-dependent average)
    expected = cooltools.expected_cis(
        clr,
        regions=[(CHROM, REGION_START, REGION_END)],
        nproc=4,
    )

    # Extract raw matrix
    raw = clr.matrix(balance=True).fetch(REGION)
    raw = np.nan_to_num(raw, nan=0.0)

    # Build expected matrix
    n = raw.shape[0]
    exp_matrix = np.zeros_like(raw)
    for i in range(n):
        for j in range(n):
            d = abs(i - j)
            if d < len(expected):
                exp_matrix[i, j] = expected.iloc[d]["balanced.avg"] if d < len(expected) else 0

    # O/E
    obs_exp = np.divide(raw, exp_matrix, where=exp_matrix > 0, out=np.zeros_like(raw))

    # P(s) curve
    max_dist = n
    ps = np.zeros(max_dist)
    counts = np.zeros(max_dist)
    for i in range(n):
        for j in range(i, n):
            d = j - i
            if d < max_dist:
                ps[d] += raw[i, j]
                counts[d] += 1
    ps = np.divide(ps, counts, where=counts > 0, out=np.zeros_like(ps))
    distances = np.arange(max_dist) * resolution

    return obs_exp, (distances, ps)


def compute_insulation_score(matrix: np.ndarray, window: int = 50) -> np.ndarray:
    """
    Compute insulation score profile for a contact matrix.

    Parameters
    ----------
    matrix : np.ndarray
        Square contact matrix.
    window : int
        Window size in bins.

    Returns
    -------
    insulation : np.ndarray
        Log2 insulation score at each position.
    """
    n = matrix.shape[0]
    insulation = np.zeros(n)

    for i in range(window, n - window):
        upstream = matrix[i - window:i, i - window:i]
        downstream = matrix[i:i + window, i:i + window]
        cross = matrix[i - window:i, i:i + window]

        mean_cross = np.nanmean(cross)
        mean_self = (np.nanmean(upstream) + np.nanmean(downstream)) / 2

        if mean_self > 0 and mean_cross > 0:
            insulation[i] = np.log2(mean_cross / mean_self)

    return insulation


def main():
    """Download or prepare Hi-C data for the Sox2 locus."""
    import argparse

    parser = argparse.ArgumentParser(description="Download/process Hi-C data")
    parser.add_argument("--output", type=str, default=None,
                        help="Output directory for data files")
    parser.add_argument("--mcool-es", type=str, default=None,
                        help="Path to mESC mcool file (if already downloaded)")
    parser.add_argument("--mcool-cn", type=str, default=None,
                        help="Path to cortical neuron mcool file")

    args = parser.parse_args()

    if args.output is None:
        args.output = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    os.makedirs(args.output, exist_ok=True)

    if args.mcool_es and args.mcool_cn:
        # Extract matrices
        logger.info("Extracting mESC matrix...")
        es_matrix = extract_region_from_mcool(args.mcool_es, RESOLUTION)
        np.save(os.path.join(args.output, "hic_mESC_Sox2.npy"), es_matrix)

        logger.info("Extracting cortical neuron matrix...")
        cn_matrix = extract_region_from_mcool(args.mcool_cn, RESOLUTION)
        np.save(os.path.join(args.output, "hic_CN_Sox2.npy"), cn_matrix)

        # Insulation scores
        ins_es = compute_insulation_score(es_matrix)
        ins_cn = compute_insulation_score(cn_matrix)
        np.save(os.path.join(args.output, "insulation_mESC.npy"), ins_es)
        np.save(os.path.join(args.output, "insulation_CN.npy"), ins_cn)

        logger.info("Data extraction complete!")
    else:
        download_from_geo(args.output)


if __name__ == "__main__":
    main()
