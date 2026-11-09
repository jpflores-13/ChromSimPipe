"""Several loci in one config: the active locus, its tiling and its folders."""
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

_PROBE = (
    "import json; from configs import parameters as p; "
    "print(json.dumps(dict(active=p.ACTIVE_LOCUS, names=p.LOCUS_NAMES, region=p.REGION, "
    "n=p.N_MONOMERS, tiling=p.TILING, ctcf=p.PATHS['ctcf'], mcool=p.PATHS['mcool'], "
    "beds=p.CTCF_BEDS)))"
)


def _config(tmp_path):
    cfg = yaml.safe_load((ROOT / "config" / "ChromSimConfig.yaml").read_text())
    cfg["loci"] = [
        dict(name="locusA", chrom="chr1", region_start=1_000_000, region_end=1_950_000),
        dict(name="locusB", chrom="chr7", region_start=5_000_000, region_end=5_500_000),
    ]
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def _load(config, locus=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CHROMSIM_")}
    env["CHROMSIM_CONFIG"] = str(config)
    if locus:
        env["CHROMSIM_LOCUS"] = locus
    res = subprocess.run([sys.executable, "-c", _PROBE], cwd=ROOT, env=env,
                         capture_output=True, text=True)
    return res


def test_first_locus_is_active_by_default(tmp_path):
    out = json.loads(_load(_config(tmp_path)).stdout)
    assert out["names"] == ["locusA", "locusB"]
    assert out["active"] == "locusA" and out["region"] == "chr1:1000000-1950000"


def test_chromsim_locus_selects_locus_and_its_tiling(tmp_path):
    out = json.loads(_load(_config(tmp_path), "locusB").stdout)
    assert out["active"] == "locusB" and out["n"] == 500
    t = out["tiling"]
    assert t["tile_size"] == 500 + 2 * t["padding"]
    assert t["chrom_size"] == t["n_tiles"] * t["tile_size"] <= 70000


def test_each_locus_has_its_own_folders_but_shares_hic(tmp_path):
    a = json.loads(_load(_config(tmp_path), "locusA").stdout)
    b = json.loads(_load(_config(tmp_path), "locusB").stdout)
    assert a["ctcf"] != b["ctcf"] and a["ctcf"].endswith(os.path.join("loci", "locusA", "ctcf"))
    assert a["mcool"] == b["mcool"]
    assert all(os.path.dirname(p) == a["ctcf"] for p in a["beds"].values())


def test_unknown_locus_is_rejected(tmp_path):
    res = _load(_config(tmp_path), "nope")
    assert res.returncode != 0 and "not one of the config's loci" in res.stderr


def _config_with(tmp_path, locus, genome=None):
    cfg = yaml.safe_load((ROOT / "config" / "ChromSimConfig.yaml").read_text())
    cfg["loci"] = [locus]
    if genome:
        cfg["genome"] = str(genome)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def test_window_at_chromosome_start_is_accepted(tmp_path):
    res = _load(_config_with(tmp_path, dict(name="tel", chrom="chr2", region_start=0,
                                            region_end=1_000_000)))
    assert res.returncode == 0 and json.loads(res.stdout)["n"] == 1000


def test_window_off_the_resolution_grid_is_rejected(tmp_path):
    res = _load(_config_with(tmp_path, dict(name="off", chrom="chr1", region_start=64_210_500,
                                            region_end=65_210_500)))
    assert res.returncode != 0 and "multiples of resolution" in res.stderr


def test_coordinates_with_units_are_rejected(tmp_path):
    res = _load(_config_with(tmp_path, dict(name="mb", chrom="chr1", region_start="64.21Mb",
                                            region_end="65.21Mb")))
    assert res.returncode != 0 and "whole numbers" in res.stderr


def test_window_past_chromosome_end_is_rejected(tmp_path):
    genome = tmp_path / "genome.fa"
    (tmp_path / "genome.fa.fai").write_text("chr1\t248956422\t6\t60\t61\n")
    res = _load(_config_with(tmp_path, dict(name="end", chrom="chr1", region_start=248_500_000,
                                            region_end=249_000_000), genome))
    assert res.returncode != 0 and "past the end of chr1" in res.stderr and "248956000" in res.stderr
