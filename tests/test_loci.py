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
