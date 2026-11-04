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


HIC2COOL_COLUMNS = ["chrom", "start", "end", "KR", "VC", "VC_SQRT"]


def test_named_normalization_is_used():
    assert exp.resolve_balance(_Cooler(HIC2COOL_COLUMNS), "KR") == "KR"


def test_raw_counts():
    assert exp.resolve_balance(_Cooler(HIC2COOL_COLUMNS), False) is False


def test_missing_weight_column_names_the_alternatives():
    with pytest.raises(ValueError, match="KR"):
        exp.resolve_balance(_Cooler(HIC2COOL_COLUMNS), True)


def test_unknown_colormap_falls_back():
    pytest.importorskip("matplotlib")
    assert exp._resolve_cmap("viridis") == "viridis"
    assert exp._resolve_cmap("not_a_cmap") == "YlOrRd"
