"""Unit tests for the flat-folder output layout in scripts/run_analysis_all.py."""
import os
import pytest


def test_build_prefix_per_replicate(tmp_path, monkeypatch):
    """Per-rep prefix combines condition + n_blocks + rep."""
    import run_analysis_all as r

    # Stub count_conformations so we don't need real .h5 files.
    monkeypatch.setattr(r, "count_conformations", lambda sim_dir: 2400)

    sim_dir = str(tmp_path / "mESC_baseline_rep0")
    os.makedirs(sim_dir, exist_ok=True)

    prefix = r._build_prefix(sim_dir)
    assert prefix == "mESC_baseline_2400blk_rep0"


def test_build_prefix_unparseable_falls_back_to_basename(tmp_path, monkeypatch):
    """When parse_condition_rep can't parse a rep number, use the basename + n_blocks."""
    import run_analysis_all as r

    monkeypatch.setattr(r, "count_conformations", lambda sim_dir: 100)
    sim_dir = str(tmp_path / "loose_dirname_no_rep")
    os.makedirs(sim_dir, exist_ok=True)

    prefix = r._build_prefix(sim_dir)
    assert prefix == "loose_dirname_no_rep_100blk"


def test_run_msd_for_pair_uses_file_prefix(tmp_path, monkeypatch):
    """When file_prefix is given, MSD JSON/NPZ files include it."""
    import msd_two_point as msd

    written = []
    monkeypatch.setattr(msd, "save_msd_json",
                        lambda path, *a, **kw: written.append(path))
    monkeypatch.setattr(msd, "save_msd_npz",
                        lambda path, *a, **kw: written.append(path))
    monkeypatch.setattr(
        msd, "compute_two_point_msd_tiled",
        lambda *a, **kw: {"lags": [1, 2, 3], "msd": [0.1, 0.2, 0.3],
                          "n_tiles_used": 1, "per_tile_msd": None},
    )
    monkeypatch.setattr(msd, "fit_msd_alpha",
                        lambda *a, **kw: {"alpha": 1.0, "K_alpha": 0.1})
    monkeypatch.setattr(msd, "fit_msd_saturation", lambda *a, **kw: {})
    monkeypatch.setattr(msd, "fit_msd_alpha_per_tile", lambda *a, **kw: {})

    msd.run_msd_for_pair(
        conformations=[None],
        pair=(0, 100),
        label="probe",
        out_dir=str(tmp_path),
        tile_size=10, pad=0, n_tiles=1,
        file_prefix="mESC_2400blk_rep0",
    )
    assert any(p.endswith("mESC_2400blk_rep0_msd_probe.json") for p in written)
    assert any(p.endswith("mESC_2400blk_rep0_msd_probe.npz") for p in written)


def test_run_msd_and_dynamics_signature():
    """_run_msd_and_dynamics must accept (common_dir, file_prefix), not out_dir."""
    import inspect
    import run_analysis_all as r

    sig = inspect.signature(r._run_msd_and_dynamics)
    assert "common_dir" in sig.parameters
    assert "file_prefix" in sig.parameters
    assert "out_dir" not in sig.parameters


def test_warn_on_legacy_unprefixed_outputs(tmp_path, caplog):
    """End-of-run scan flags un-prefixed legacy filenames in common_dir."""
    import logging
    import run_analysis_all as r

    # Simulate a botched run that left these around.
    (tmp_path / "sim_contact_map.npy").write_bytes(b"")
    (tmp_path / "msd_overlay.png").write_bytes(b"")
    (tmp_path / "mESC_2400blk_rep0_contact_map.npy").write_bytes(b"")  # OK

    with caplog.at_level(logging.WARNING, logger=r.logger.name):
        r._warn_legacy_outputs(str(tmp_path))

    msgs = [rec.getMessage() for rec in caplog.records]
    assert any("sim_contact_map.npy" in m for m in msgs)
    assert any("msd_overlay.png" in m for m in msgs)
    assert not any("mESC_2400blk_rep0_contact_map.npy" in m for m in msgs)


def test_parse_condition_rep_strips_merged_prefix():
    """A 'merged_' prefix on the dir name is dropped before condition parsing,
    so output filenames stay condition-only (no 'merged_' bleed-through)."""
    import run_analysis_all as r

    legacy = r.parse_condition_rep("mESC_ctcf-mESC_rep0")
    new = r.parse_condition_rep("merged_mESC_ctcf-mESC_rep0")
    assert legacy == ("mESC_ctcf-mESC", 0)
    assert new == ("mESC_ctcf-mESC", 0)


def test_build_prefix_for_merged_dir(tmp_path, monkeypatch):
    """_build_prefix on a merged_-prefixed sim dir yields the same stem as legacy."""
    import run_analysis_all as r
    monkeypatch.setattr(r, "count_conformations", lambda sim_dir: 2400)

    sim_dir = str(tmp_path / "merged_mESC_baseline_rep0")
    os.makedirs(sim_dir, exist_ok=True)
    assert r._build_prefix(sim_dir) == "mESC_baseline_2400blk_rep0"


def test_migrate_one_renames_legacy_files(tmp_path, monkeypatch):
    """migrate_legacy_analysis moves legacy <sim_dir>/analysis/ files
    into the flat folder with the same renaming the live pipeline applies.
    """
    import migrate_legacy_analysis as m
    import run_analysis_all as r

    monkeypatch.setattr(r, "count_conformations", lambda sim_dir: 2400)

    sim_dir = tmp_path / "polychrom_3d" / "mESC_baseline_rep0"
    legacy = sim_dir / "analysis"
    legacy.mkdir(parents=True)
    (legacy / "sim_contact_map.npy").write_bytes(b"x")
    (legacy / "msd_overlay.png").write_bytes(b"y")
    (legacy / "comparison_metrics.json").write_text("{}")
    (legacy / "msd_probeA.json").write_text("{}")  # per-pair MSD artefact

    out = tmp_path / "results" / "analysis"

    n = m.migrate_one(str(sim_dir), str(out),
                      copy=False, dry_run=False, delete_empty=True)
    assert n == 4
    assert (out / "mESC_baseline_2400blk_rep0_contact_map.npy").exists()
    assert (out / "mESC_baseline_2400blk_rep0_msd_overlay.png").exists()
    assert (out / "mESC_baseline_2400blk_rep0_metrics.json").exists()
    assert (out / "mESC_baseline_2400blk_rep0_msd_probeA.json").exists()
    # Legacy folder removed because it was emptied.
    assert not legacy.exists()


def test_migrate_one_dry_run_does_not_touch_disk(tmp_path, monkeypatch):
    """Dry-run logs intent but leaves the source intact and never creates dst."""
    import migrate_legacy_analysis as m
    import run_analysis_all as r

    monkeypatch.setattr(r, "count_conformations", lambda sim_dir: 100)

    sim_dir = tmp_path / "polychrom_3d" / "CN_baseline_rep1"
    legacy = sim_dir / "analysis"
    legacy.mkdir(parents=True)
    (legacy / "sim_contact_map.npy").write_bytes(b"x")
    out = tmp_path / "results" / "analysis"

    n = m.migrate_one(str(sim_dir), str(out),
                      copy=False, dry_run=True, delete_empty=True)
    assert n == 1
    assert (legacy / "sim_contact_map.npy").exists()
    assert not out.exists()


def test_plot_panels_discover_maps(tmp_path):
    """plot_contact_map_panels.discover_maps groups per-rep and pooled files."""
    import plot_contact_map_panels as p

    d = tmp_path / "analysis"
    d.mkdir()
    # Per-replicate files for two conditions.
    (d / "mESC_baseline_2400blk_rep0_contact_map.npy").write_bytes(b"")
    (d / "mESC_baseline_2400blk_rep1_contact_map.npy").write_bytes(b"")
    (d / "CN_short_residency_3000blk_rep0_contact_map.npy").write_bytes(b"")
    # Pooled file for one condition.
    (d / "mESC_baseline_4800blk_pooled_contact_map.npy").write_bytes(b"")
    # Decoy files that must NOT be classified as contact maps.
    (d / "mESC_baseline_2400blk_rep0_msd_overlay.png").write_bytes(b"")
    (d / "ps_overlay_all_conditions.png").write_bytes(b"")

    per_rep, pooled = p.discover_maps(str(d))
    assert set(per_rep) == {"mESC_baseline", "CN_short_residency"}
    assert [r for r, _ in per_rep["mESC_baseline"]] == [0, 1]
    assert [r for r, _ in per_rep["CN_short_residency"]] == [0]
    assert set(pooled) == {"mESC_baseline"}
    assert pooled["mESC_baseline"][0] == 4800
