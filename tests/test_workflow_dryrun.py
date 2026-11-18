"""The Snakemake workflow: dry runs and the settings that trigger reruns.

The full dry run needs the real inputs the test config names (genome, CTCF
peaks, Hi-C, linked into data/ per data/README.md); it is skipped when any of
them is absent. Everything else here runs from the config and samplesheet alone.
"""
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "workflows" / "ChromSimPipe.snakefile"
TEST_CONFIG = ROOT / "config" / "ChromSimConfig.test.yaml"
FULL_CONFIG = ROOT / "config" / "ChromSimConfig.yaml"


def _clean_env(**extra):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CHROMSIM_")}
    env.update(extra)
    return env


def _probe(config, code, overrides=None):
    """Run `code` against configs.parameters loaded from `config`; return its JSON output."""
    env = _clean_env(CHROMSIM_CONFIG=str(config))
    if overrides is not None:
        env["CHROMSIM_CONFIG_OVERRIDES"] = json.dumps(overrides)
    r = subprocess.run([sys.executable, "-c",
                        "import json; from configs import parameters as p; " + code],
                       cwd=ROOT, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return json.loads(r.stdout)


def _job_params(config, stage, overrides=None):
    """{condition: job_params_json} at the first locus, for one stage."""
    return _probe(config, (
        f"print(json.dumps({{c['name']: p.job_params_json({stage!r}, p.LOCUS_NAMES[0], c['name']) "
        f"for c in p.SIMULATION_CONDITIONS}}))"), overrides)


def _edited_config(tmp_path, edit):
    cfg = yaml.safe_load(FULL_CONFIG.read_text())
    edit(cfg)
    path = tmp_path / "edited.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def _snakemake():
    exe = shutil.which("snakemake")
    if exe is None:
        pytest.skip("snakemake is not installed (pip install -r config/requirements.txt)")
    return exe


def _dry_run(tmp_path, *args):
    """snakemake -n on the test config, with every output under tmp_path."""
    r = subprocess.run(
        [_snakemake(), "-s", str(SNAKEFILE), "--configfile", str(TEST_CONFIG),
         "-n", "--cores", "1", *args],
        cwd=ROOT, env=_clean_env(CHROMSIM_OUTPUT_ROOT=str(tmp_path)),
        capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def _job_stats(log):
    """{rule: count} from the 'Job stats' table of a dry run."""
    table = log.split("Job stats:", 1)[1]
    counts = {}
    for line in table.splitlines()[1:]:
        m = re.match(r"^(\S+)\s+(\d+)$", line.strip())
        if m:
            counts[m.group(1)] = int(m.group(2))
        elif counts and not line.strip():
            break
    counts.pop("total", None)
    return counts


# ---------------------------------------------------------------------------
# Settings reach configs.parameters and the jobs' params
# ---------------------------------------------------------------------------

def test_config_overrides_are_merged():
    """`--config` values (passed on as $CHROMSIM_CONFIG_OVERRIDES) override the
    file, nested keys merge, and untouched keys keep their file values."""
    base = _probe(TEST_CONFIG, "print(json.dumps([p.N_SHARDS, p.POLYSIM]))")
    n_shards, polysim = _probe(TEST_CONFIG, "print(json.dumps([p.N_SHARDS, p.POLYSIM]))",
                               overrides={"n_shards": base[0] + 3,
                                          "polysim": {"vlef": base[1]["vlef"] * 2}})
    assert n_shards == base[0] + 3
    assert polysim["vlef"] == base[1]["vlef"] * 2
    assert {k: v for k, v in polysim.items() if k != "vlef"} == \
           {k: v for k, v in base[1].items() if k != "vlef"}


@pytest.mark.parametrize("stage", ["simulate", "merge", "analyze"])
def test_editing_one_condition_changes_only_its_params(tmp_path, stage):
    before = _job_params(FULL_CONFIG, stage)
    if len(before) < 2:
        pytest.skip("needs a config with at least two conditions")
    first = next(iter(before))

    def edit(cfg):
        cfg["conditions"][0]["lifetime"] = int(cfg["conditions"][0]["lifetime"]) * 2

    after = _job_params(_edited_config(tmp_path, edit), stage)
    assert after[first] != before[first]
    assert {k: v for k, v in after.items() if k != first} == \
           {k: v for k, v in before.items() if k != first}


def test_n_shards_override_changes_every_simulation_job():
    n_shards = _probe(FULL_CONFIG, "print(json.dumps(p.N_SHARDS))")
    before = _job_params(FULL_CONFIG, "simulate")
    after = _job_params(FULL_CONFIG, "simulate", overrides={"n_shards": n_shards + 2})
    assert all(after[c] != before[c] for c in before)


def test_analysis_settings_change_only_analysis_params(tmp_path):
    def edit(cfg):
        cfg["analysis"]["apa"]["pad"] = int(cfg["analysis"]["apa"]["pad"]) + 1

    path = _edited_config(tmp_path, edit)
    assert _job_params(path, "simulate") == _job_params(FULL_CONFIG, "simulate")
    assert _job_params(path, "merge") == _job_params(FULL_CONFIG, "merge")
    after, before = _job_params(path, "analyze"), _job_params(FULL_CONFIG, "analyze")
    assert all(after[c] != before[c] for c in before)


def test_descriptions_and_labels_rerun_nothing(tmp_path):
    def edit(cfg):
        for c in cfg["conditions"]:
            c["description"] = "edited " + str(c.get("description", ""))
        for loc in cfg["loci"]:
            loc["label"] = "edited " + str(loc.get("label", ""))

    path = _edited_config(tmp_path, edit)
    for stage in ("simulate", "merge", "analyze"):
        assert _job_params(path, stage) == _job_params(FULL_CONFIG, stage)


# ---------------------------------------------------------------------------
# Snakemake dry runs
# ---------------------------------------------------------------------------

_LOCALRULES = """
import json, sys
from pathlib import Path
from snakemake.api import SnakemakeApi
from snakemake.settings.types import ConfigSettings, OutputSettings, ResourceSettings
with SnakemakeApi(OutputSettings(quiet={"all"})) as api:
    wf = api.workflow(snakefile=Path(sys.argv[1]), workdir=Path(sys.argv[2]),
                      config_settings=ConfigSettings(configfiles=[Path(sys.argv[3])]),
                      resource_settings=ResourceSettings())
    w = wf._workflow
    print(json.dumps(dict(rules=sorted(r.name for r in w.rules),
                          local=sorted(r.name for r in w.rules if r.name in w._localrules))))
"""


def test_getmotif_is_a_localrule(tmp_path):
    """getMotif downloads with wget, so it must run on the node running
    Snakemake (compute nodes may have no outside network). Needs no inputs."""
    pytest.importorskip("snakemake")
    motif = _probe(TEST_CONFIG, "print(json.dumps(p._cfg['ctcf_motif']))")
    if motif.get("meme"):
        pytest.skip("the test config uses a local MEME file: no getMotif rule")
    r = subprocess.run([sys.executable, "-c", _LOCALRULES, str(SNAKEFILE), str(ROOT), str(TEST_CONFIG)],
                       cwd=ROOT, env=_clean_env(CHROMSIM_OUTPUT_ROOT=str(tmp_path)),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rules = json.loads(r.stdout.strip().splitlines()[-1])
    assert "getMotif" in rules["rules"]
    assert rules["local"] == ["getMotif"]


def test_snakefile_dry_run_without_inputs(tmp_path):
    """The Snakefile parses and plans the rules that need no input files."""
    motif = _probe(TEST_CONFIG, "print(json.dumps(p._cfg['ctcf_motif']))")
    out = Path(_probe(TEST_CONFIG, "print(json.dumps(p.OUTPUT_DIR))")).relative_to(ROOT)
    targets = [str(tmp_path / out / "done" / "preparePolysim.done")]
    if not motif.get("meme"):
        targets.append(str(tmp_path / out / "motifs" / f"{motif['jaspar_id']}.meme"))
    code, log = _dry_run(tmp_path, *targets)
    assert code == 0, log
    assert _job_stats(log) == ({"preparePolysim": 1} if motif.get("meme")
                               else {"preparePolysim": 1, "getMotif": 1})


def _required_inputs(config):
    """Every input file outside the workflow's outputs that a full run reads."""
    return _probe(config, (
        "conds = p.SIMULATION_CONDITIONS; "
        "files = [p.GENOME] + [p.CTCF_PEAKS[t] for t in p.CTCF_TYPES] "
        "+ [p.HIC_FILES[c['compare_hic']] for c in conds if c['compare_hic']] "
        "+ [f for c in conds for f in p.condition_peak_files(c).values()] "
        "+ ([p._repo_path(p._cfg['ctcf_motif']['meme'])] if p._cfg['ctcf_motif'].get('meme') else []); "
        "print(json.dumps(sorted(set(files))))"))


def _expected_counts(config):
    """Job count per rule, derived from the config."""
    c = _probe(config, (
        "print(json.dumps(dict("
        "loci=len(p.LOCUS_NAMES), conds=len(p.SIMULATION_CONDITIONS), "
        "reps=p.N_REPLICATES, shards=p.N_SHARDS, ctcf_types=len(p.CTCF_TYPES), "
        "hic=len({c['compare_hic'] for c in p.SIMULATION_CONDITIONS if c['compare_hic']}), "
        "local_meme=bool(p._cfg['ctcf_motif'].get('meme')))))"))
    counts = dict(
        all              = 1,
        convertHic       = c["hic"],
        getMotif         = 0 if c["local_meme"] else 1,
        preparePolysim   = 1,
        orientCTCF       = c["loci"] * c["ctcf_types"],
        validateCTCF     = c["loci"],
        simulateShard    = c["loci"] * c["conds"] * c["reps"] * c["shards"],
        mergeShards      = c["loci"] * c["conds"] * c["reps"],
        analyzeCondition = c["loci"] * c["conds"],
        benchmarkSummary = 1,
    )
    return {k: v for k, v in counts.items() if v}


def test_dry_run_job_counts(tmp_path):
    missing = [f for f in _required_inputs(TEST_CONFIG) if not os.path.exists(f)]
    if missing:
        pytest.skip(f"real inputs absent (see data/README.md): {', '.join(missing)}")
    code, log = _dry_run(tmp_path)
    assert code == 0, log
    assert _job_stats(log) == _expected_counts(TEST_CONFIG)


def test_dry_run_job_counts_follow_config_overrides(tmp_path):
    """--config n_shards / n_replicates change how many jobs Snakemake plans."""
    missing = [f for f in _required_inputs(TEST_CONFIG) if not os.path.exists(f)]
    if missing:
        pytest.skip(f"real inputs absent (see data/README.md): {', '.join(missing)}")
    base = _expected_counts(TEST_CONFIG)
    code, log = _dry_run(tmp_path, "--config", "n_shards=2", "n_replicates=2")
    assert code == 0, log
    stats = _job_stats(log)
    per = base["analyzeCondition"]
    assert stats["simulateShard"] == per * 2 * 2
    assert stats["mergeShards"] == per * 2
