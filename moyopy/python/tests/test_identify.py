from __future__ import annotations

import pytest

from moyopy import PointGroup, SpaceGroup, operations_from_number


@pytest.mark.parametrize(
    "number,arithmetic_number",
    [
        (158, 45),  # P3c1, 3m1P
        (159, 46),  # P31c, 31mP
    ],
)
def test_identify_point_group(number: int, arithmetic_number: int):
    operations = operations_from_number(number)
    point_group = PointGroup(operations.rotations)
    assert point_group.arithmetic_number == arithmetic_number


@pytest.mark.parametrize(
    "number",
    [
        158,  # P3c1
        159,  # P31c
    ],
)
def test_identify_space_group(number: int):
    operations = operations_from_number(number)
    space_group = SpaceGroup(
        prim_rotations=operations.rotations, prim_translations=operations.translations
    )
    assert space_group.number == number
