"""FIMO hits are matched back to peaks whatever form FIMO reports them in.

FIMO_559 is verbatim `fimo --text` output from MEME 5.5.9 on Longleaf, run on
three STRS CTCF peak windows made by `bedtools getfasta -name` (v2.31.1).
FIMO reported genome coordinates on `chr1`, not the window names.
"""
import extract_ctcf_sites as ext

FIMO_559 = "\n".join([
    "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence",
    "MA0139.1\tCTCF\tchr1\t64248197\t64248215\t-\t11.3443\t2.04e-05\t\tGGGACAGGAGAGGGAGCAA",
    "MA0139.1\tCTCF\tchr1\t64469999\t64470017\t-\t3.70492\t0.000467\t\tGGGACGCCGGGGGCCGGGG",
    "MA0139.1\tCTCF\tchr1\t64470152\t64470170\t+\t3.01639\t0.000593\t\tGCGCAGCGTGGTGGCACCA",
])

SUMMITS = [64248123, 64269369, 64470086]
FLANK = 100
NAME_TO_COORD = {f"chr1_{s}": ("chr1", s) for s in SUMMITS}
WINDOWS = [("chr1", s - FLANK, s + FLANK, ("chr1", s)) for s in SUMMITS]


def test_genome_coordinate_output_is_matched():
    hits = ext._best_hits_from_fimo(FIMO_559, NAME_TO_COORD, WINDOWS)
    assert hits[("chr1", 64248123)] == ("-", 2.04e-05)
    # two hits in the third window: the lower p-value wins
    assert hits[("chr1", 64470086)] == ("-", 0.000467)
    # no hit at p < 1e-3 in the second window
    assert ("chr1", 64269369) not in hits


def test_window_relative_output_is_matched():
    """Same hit reported against a `chr:start-end` header (window-relative)."""
    start0 = 64248123 - FLANK
    relative = FIMO_559.splitlines()[0] + "\n" + "\t".join([
        "MA0139.1", "CTCF", f"chr1:{start0}-{start0 + 2 * FLANK}",
        str(64248197 - start0), str(64248215 - start0), "-", "11.3443", "2.04e-05", "", "X",
    ])
    hits = ext._best_hits_from_fimo(relative, NAME_TO_COORD, WINDOWS)
    assert hits == {("chr1", 64248123): ("-", 2.04e-05)}
