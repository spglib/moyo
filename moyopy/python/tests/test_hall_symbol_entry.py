from __future__ import annotations

from moyopy import HallSymbolEntry


def test_hall_symbol_entry():
    hs = HallSymbolEntry(hall_number=528)  # No. 228, origin choice 2
    assert hs.hall_number == 528
    assert hs.number == 228
    assert hs.setting == "2"
    assert hs.hm_short == "F d -3 c"
