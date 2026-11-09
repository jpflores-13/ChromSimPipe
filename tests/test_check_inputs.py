"""Helpers of scripts/check_inputs.py (no input files are read)."""
import check_inputs as ci


def test_other_naming_style():
    assert ci.other_style("chr7") == "7"
    assert ci.other_style("7") == "chr7"


def test_hint_names_the_other_style_when_present():
    assert "'7'" in ci.name_hint("chr7", {"7", "8"})
    assert ci.name_hint("chr7", {"chr8"}) == ""


def test_genome_sizes_come_from_the_index(tmp_path):
    (tmp_path / "g.fa.fai").write_text("chr1\t248956422\t6\t60\t61\n")
    assert ci.genome_sizes(str(tmp_path / "g.fa")) == {"chr1": 248956422}
