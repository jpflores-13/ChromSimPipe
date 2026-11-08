"""polysim wiring: config parsing and the SimParams arguments built for a shard.

These check ChromSimPipe's side of the integration (orientation mapping, per-site
stall probabilities, schedule validation). They do not import polysim itself,
which needs OpenMM / polychrom.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

_PROBE = r"""
import json, numpy as np
from scripts import run_simulation_shard as rs
from configs.parameters import get_condition, POLYSIM
cond = get_condition("cond1_baseline")
kw = rs.polysim_kwargs(cond["params"], np.array([10, 20, 30]), np.array([1, -1, 1]),
                       n_frames=FRAMES, gpu=0, shard_dir="/tmp/x")
print(json.dumps(dict(polysim=POLYSIM, kw={k: v for k, v in kw.items()
     if k not in ("ctcf_left", "ctcf_right")}, left=kw["ctcf_left"], right=kw["ctcf_right"])))
"""


def _run(frames=2000):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CHROMSIM_")}
    env["CHROMSIM_CONFIG"] = str(ROOT / "config" / "ChromSimConfig.test.yaml")
    res = subprocess.run([sys.executable, "-c", _PROBE.replace("FRAMES", str(frames))],
                         cwd=ROOT, env=env, capture_output=True, text=True)
    return res


def test_polysim_settings_parse():
    res = _run()
    assert res.returncode == 0, res.stderr
    out = json.loads(res.stdout)
    assert out["polysim"]["vlef"] == 0.0025
    assert out["polysim"]["life_boost_stalled"] == 4.0


def test_orientation_and_stall_probabilities():
    out = json.loads(_run().stdout)
    # +1 sites stop leftward-moving legs (polysim ctcf_left), -1 the rightward ones
    assert out["left"] == {"10": 0.5, "30": 0.5}
    assert out["right"] == {"20": 0.5}
    kw = out["kw"]
    assert kw["life"] == 75000 and kw["sep"] == 240 and kw["stall"] is None
    assert kw["numsave"] == 2000 and kw["saveevery"] == 50
    assert kw["lifebooststalled"] == 4.0 and kw["vlef"] == 0.0025


def test_schedule_mismatch_is_reported_up_front():
    res = _run(frames=1990)  # 1990 x 50 blocks is not a multiple of 1000
    assert res.returncode != 0
    assert "multiple of 20 frames per shard" in res.stderr
