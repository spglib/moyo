from __future__ import annotations

from moyopy import ArithmeticCrystalClass


def test_arithmetic_crystal_class():
    a1 = ArithmeticCrystalClass(45)  # 3m1P
    assert a1.arithmetic_number == 45
    assert a1.symbol == "3m1P"
    assert a1.geometric_crystal_class == "3m"
    assert a1.bravais_class == "hP"
