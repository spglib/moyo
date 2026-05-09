from __future__ import annotations

from moyopy import LayerHallSymbolEntry


def test_layer_hall_symbol_entry_first():
    entry = LayerHallSymbolEntry(hall_number=1)
    assert entry.hall_number == 1
    assert entry.number == 1
    assert entry.arithmetic_number == 1
    assert entry.setting == ""
    assert entry.hall_symbol == "p 1"
    assert entry.hm_short == "p 1"
    assert entry.hm_full == "p 1"
    assert entry.centering.order == 1


def test_layer_hall_symbol_entry_last():
    entry = LayerHallSymbolEntry(hall_number=116)
    assert entry.hall_number == 116
    # LG 80 / p6/mmm
    assert entry.number == 80
    assert entry.arithmetic_number == 43
    assert entry.hm_short == "p 6/m m m"
