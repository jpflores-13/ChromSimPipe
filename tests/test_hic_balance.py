"""Choosing the Hi-C normalization column, and colormap fallback."""
import pytest

import experimental_compare as exp


class _Bins:
    def __init__(self, columns):
        self.columns = columns


class _Cooler:
    """Stand-in exposing only what resolve_balance reads (no contact data)."""
    uri = "test.cool"

    def __init__(self, columns):
        self._columns = columns

    def bins(self):
        return _Bins(self._columns)


## Bins columns hictk convert writes for a .hic: one per normalization stored.
V8_COLUMNS = ["chrom", "start", "end", "KR", "VC", "VC_SQRT"]          # Juicer 1 (.hic v8)
V9_COLUMNS = ["chrom", "start", "end", "SCALE", "VC", "VC_SQRT"]       # Juicer 2 (.hic v9)


def test_named_normalization_is_used():
    assert exp.resolve_balance(_Cooler(V8_COLUMNS), "KR") == "KR"


def test_raw_counts():
    assert exp.resolve_balance(_Cooler(V8_COLUMNS), False) is False


def test_missing_weight_column_names_the_alternatives():
    with pytest.raises(ValueError, match="KR"):
        exp.resolve_balance(_Cooler(V8_COLUMNS), True)


def test_auto_prefers_kr_from_converted_hic():
    assert exp.resolve_balance(_Cooler(V8_COLUMNS), "auto") == "KR"


def test_auto_uses_weight_from_balanced_mcool():
    assert exp.resolve_balance(_Cooler(["chrom", "start", "end", "weight"]), "auto") == "weight"


def test_auto_falls_back_to_scale_without_kr():
    assert exp.resolve_balance(_Cooler(V9_COLUMNS), "auto") == "SCALE"


def test_auto_prefers_kr_over_scale():
    assert exp.resolve_balance(_Cooler(V9_COLUMNS + ["KR"]), "auto") == "KR"


def test_scale_by_name():
    assert exp.resolve_balance(_Cooler(V9_COLUMNS), "SCALE") == "SCALE"


def test_pick_balance_on_names_from_a_hic_header():
    assert exp.pick_balance(["SCALE", "VC", "VC_SQRT"], "auto", "x.hic") == "SCALE"
    with pytest.raises(ValueError, match="x.hic"):
        exp.pick_balance(["VC"], "KR", "x.hic")


def test_auto_without_kr_scale_or_weight_raises():
    with pytest.raises(ValueError, match="auto"):
        exp.resolve_balance(_Cooler(["chrom", "start", "end", "VC"]), "auto")


def test_unknown_colormap_falls_back():
    pytest.importorskip("matplotlib")
    assert exp._resolve_cmap("viridis") == "viridis"
    assert exp._resolve_cmap("not_a_cmap") == "YlOrRd"
