"""configs/parameters.py builds everything from the config + samplesheet."""
import os

import pytest


def test_conditions_resolve_to_samplesheet_rows():
    from configs import parameters as p

    for cond in p.SIMULATION_CONDITIONS:
        assert cond["ctcf_type"] in p.CTCF_PEAKS
        if cond["compare_hic"]:
            assert cond["compare_hic"] in p.HIC_FILES


def test_outputs_live_under_output_dir():
    from configs import parameters as p

    paths = [*p.PATHS.values(), *p.CTCF_BEDS.values(), *p.MCOOL_FILES.values()]
    for path in paths:
        assert os.path.commonpath([path, p.OUTPUT_DIR]) == p.OUTPUT_DIR


def test_missing_ctcf_bed_raises_instead_of_placeholder(monkeypatch):
    from configs import parameters as p

    ctcf_type = p.CTCF_TYPES[0]
    monkeypatch.delitem(p.CTCF_SITES_BY_TYPE, ctcf_type, raising=False)
    with pytest.raises(FileNotFoundError):
        p.get_ctcf_arrays(ctcf_type)


def test_tiling_matches_locus():
    from configs import parameters as p

    assert p.TILING["tile_size"] == p.N_MONOMERS + 2 * p.TILING["padding"]


def _conditions_config(tmp_path, conditions, calibration=None):
    import yaml
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config" / "ChromSimConfig.yaml").read_text())
    cfg["conditions"] = conditions
    if calibration:
        cfg["analysis"]["calibration"]["msd"].update(calibration)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


def _probe(config, code):
    import os, subprocess, sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    env = {k: v for k, v in os.environ.items() if not k.startswith("CHROMSIM_")}
    env["CHROMSIM_CONFIG"] = str(config)
    return subprocess.run([sys.executable, "-c", "import json; from configs import parameters as p; " + code],
                          cwd=root, env=env, capture_output=True, text=True)


BASE = dict(ctcf_type="condition1", compare_hic="condition1")


def test_new_condition_options_parse(tmp_path):
    import json
    cfg = _conditions_config(tmp_path, [
        dict(name="no-cohesin", separation="none", **BASE),
        dict(name="by-signal", lifetime=75000, separation=240,
             ctcf_capture={"min": 0.2, "max": 0.9}, attraction_energy=0.1, **BASE),
    ])
    res = _probe(cfg, "print(json.dumps([c['params'] for c in p.SIMULATION_CONDITIONS]))")
    assert res.returncode == 0, res.stderr
    a, b = json.loads(res.stdout)
    assert a["separation"] is None
    assert b["ctcf_capture"] == {"min": 0.2, "max": 0.9, "by": "signal"} and b["attraction_energy"] == 0.1


def test_sticky_needs_an_existing_samplesheet_column(tmp_path):
    cfg = _conditions_config(tmp_path, [
        dict(name="yy1", lifetime=75000, separation=240, ctcf_capture=0.5,
             sticky=[{"factor": "YY1", "peaks": "YY1_Peaks", "energy": 2.0}], **BASE)])
    res = _probe(cfg, "")
    assert res.returncode != 0 and "no 'YY1_Peaks' column" in res.stderr


def test_reference_condition_must_exist(tmp_path):
    cfg = _conditions_config(tmp_path, [dict(name="a", lifetime=75000, separation=240,
                                             ctcf_capture=0.5, **BASE)],
                             calibration={"reference_condition": "nope"})
    res = _probe(cfg, "")
    assert res.returncode != 0 and "reference_condition" in res.stderr
