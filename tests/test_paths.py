"""Path resolution in configs/parameters.py: samplesheet and genome paths
relative to the repo root, outputDir relative to $CHROMSIM_OUTPUT_ROOT."""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "ChromSimConfig.yaml"

_PROBE = (
    "import json; from configs import parameters as p; "
    "print(json.dumps(dict(out=p.OUTPUT_DIR, root=p.OUTPUT_ROOT, genome=p.GENOME, "
    "peaks=p.CTCF_PEAKS, hic=p.HIC_FILES)))"
)


def _load(config, **env):
    """Import parameters in a fresh interpreter (it reads env vars at import)."""
    full = {k: v for k, v in os.environ.items() if not k.startswith("CHROMSIM_")}
    full.update(CHROMSIM_CONFIG=str(config), **env)
    res = subprocess.run([sys.executable, "-c", _PROBE], cwd=ROOT, env=full,
                         capture_output=True, text=True, check=True)
    return json.loads(res.stdout)


def _config_with(tmp_path, **overrides):
    import yaml

    cfg = yaml.safe_load(CONFIG.read_text())
    cfg.update(overrides)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def test_samplesheet_inputs_resolve_under_repo_data(tmp_path):
    p = _load(_config_with(tmp_path))
    for path in [p["genome"], *p["peaks"].values(), *p["hic"].values()]:
        assert os.path.isabs(path)
        assert Path(path).parent == ROOT / "data"


def test_output_dir_defaults_to_repo_root(tmp_path):
    p = _load(_config_with(tmp_path, outputDir="output"))
    assert p["out"] == str(ROOT / "output")


def test_output_root_moves_relative_output_dir(tmp_path):
    scratch = tmp_path / "scratch"
    p = _load(_config_with(tmp_path, outputDir="output_test"),
              CHROMSIM_OUTPUT_ROOT=str(scratch))
    assert p["root"] == str(scratch)
    assert p["out"] == str(scratch / "output_test")


def test_absolute_output_dir_ignores_output_root(tmp_path):
    absolute = tmp_path / "elsewhere"
    p = _load(_config_with(tmp_path, outputDir=str(absolute)),
              CHROMSIM_OUTPUT_ROOT=str(tmp_path / "scratch"))
    assert p["out"] == str(absolute)
