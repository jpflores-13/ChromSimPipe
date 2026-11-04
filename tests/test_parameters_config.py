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
