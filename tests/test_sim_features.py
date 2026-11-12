"""Per-condition options turned into polysim inputs (scripts/sim_features.py)."""
import numpy as np
import pytest

import sim_features as feat
import sweep_1d

TILING = dict(n_tiles=3, tile_size=100, padding=10)


def test_tile_repeats_positions_in_every_copy():
    assert feat.tile([0, 5], TILING).tolist() == [10, 15, 110, 115, 210, 215]


def test_capture_by_signal_spans_min_to_max_by_rank():
    p = feat.capture_by_signal([30.0, 10.0, 20.0], 0.2, 0.8)
    assert np.allclose(p, [0.8, 0.2, 0.5])


def test_capture_by_signal_ties_share_a_rank():
    p = feat.capture_by_signal([5.0, 5.0, 9.0], 0.0, 1.0)
    assert p[0] == p[1] and p[2] == 1.0


def test_stall_dicts_follow_orientation_and_keep_the_strongest():
    left, right = feat.stall_dicts([4, 4, 9], [1, 1, -1], [0.3, 0.7, 0.5])
    assert left == {4: 0.7} and right == {9: 0.5}


def test_sticky_types_combine_factors_and_cross_energies():
    types, matrix, labels = feat.sticky_types(
        6, {"A": [1, 2], "B": [2, 4]}, {"A": 2.0, "B": 3.0},
        [{"factors": ["A", "B"], "energy": 1.0}])
    assert labels == ["none", "A", "A+B", "B"]
    assert types.tolist() == [0, 1, 2, 0, 3, 0]
    assert matrix[1, 1] == 2.0 and matrix[3, 3] == 3.0 and matrix[1, 3] == 1.0
    assert matrix[2, 2] == 3.0 and np.allclose(matrix, matrix.T) and not matrix[0].any()


def test_read_peaks_by_summit_and_by_span(tmp_path):
    bed = tmp_path / "peaks.narrowPeak"
    bed.write_text("chr1\t1000\t3000\tp1\t50\t.\t7.5\t-1\t-1\t500\n"
                   "chr2\t1000\t3000\tp2\t50\t.\t9.0\t-1\t-1\t500\n")
    m, sig = feat.read_peaks(str(bed), "chr1", 0, 10000, 1000, "summit")
    assert m.tolist() == [1] and sig[1] == 7.5
    m, _ = feat.read_peaks(str(bed), "chr1", 0, 10000, 1000, "peak")
    assert m.tolist() == [1, 2]


def test_lef_profiles_change_loading_and_stepping_only_inside_the_block():
    pytest.importorskip("polysim")
    from polysim import extrusion
    n = 50
    stall = np.zeros(n)
    with feat.lef_profiles(([10], 5.0), ([20], 4.0)):
        a = extrusion.build_lef_arrays(n, lifetime=100, vlef=0.01, stall_left=stall, stall_right=stall)
    b = extrusion.build_lef_arrays(n, lifetime=100, vlef=0.01, stall_left=stall, stall_right=stall)
    assert a["birth"][10] == 5 * b["birth"][10] and a["birth"][11] == b["birth"][11]
    assert np.isclose(1 - a["pause"][20], (1 - b["pause"][20]) / 4)


def test_sweep_scores_a_cohesin_holding_both_anchors():
    a_all, b_all = np.array([[2, 12]]), np.array([[6, 16]])   # one loop, two copies
    looped, size, cover = sweep_1d.score_frame([2, 13], [6, 15], a_all, b_all, 0, 20)
    assert looped.tolist() == [0.5] and size == 3.0 and np.isclose(cover, 6 / 20)


def test_stall_dicts_default_keeps_highest_per_monomer():
    left, right = feat.stall_dicts([5, 5, 9], [+1, +1, -1], [0.3, 0.6, 0.4])
    assert left == {5: 0.6} and right == {9: 0.4}


def test_stall_dicts_independent_combines_distinct_peaks():
    # two peaks on one monomer: 1 - (1 - 0.5)(1 - 0.5) = 0.75 (Harvey's 8.3 with n = 2)
    left, _ = feat.stall_dicts([5, 5], [+1, +1], [0.5, 0.5], groups=["pA", "pB"], combine="independent")
    assert left[5] == pytest.approx(0.75)
    # three sites at 0.5 reproduce the paper's 3xCBS value
    left, _ = feat.stall_dicts([1, 1, 1], [+1] * 3, [0.5] * 3, groups="abc", combine="independent")
    assert left[1] == pytest.approx(0.875)


def test_stall_dicts_independent_counts_one_peak_once():
    left, _ = feat.stall_dicts([5, 5], [+1, +1], [0.4, 0.5], groups=["pA", "pA"], combine="independent")
    assert left[5] == pytest.approx(0.5)


def test_stall_dicts_keeps_orientations_apart():
    left, right = feat.stall_dicts([5, 5], [+1, -1], [0.5, 0.5], groups=["a", "b"], combine="independent")
    assert left == {5: 0.5} and right == {5: 0.5}
