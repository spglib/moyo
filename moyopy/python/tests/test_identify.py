from __future__ import annotations

import pytest

from moyopy import (
    MagneticSpaceGroup,
    PointGroup,
    SpaceGroup,
    magnetic_operations_from_uni_number,
    operations_from_number,
)


@pytest.mark.parametrize(
    "number,arithmetic_number",
    [
        (158, 45),  # P3c1, 3m1P
        (159, 46),  # P31c, 31mP
    ],
)
def test_identify_point_group(number: int, arithmetic_number: int):
    operations = operations_from_number(number, primitive=True)
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
    operations = operations_from_number(number, primitive=True)
    space_group = SpaceGroup(
        prim_rotations=operations.rotations, prim_translations=operations.translations
    )
    assert space_group.number == number


@pytest.mark.parametrize(
    "uni_number",
    [
        1242,  # R31'_c
    ],
)
def test_identify_magnetic_space_group(uni_number: int):
    magnetic_operations = magnetic_operations_from_uni_number(uni_number, primitive=True)
    magnetic_space_group = MagneticSpaceGroup(
        prim_rotations=magnetic_operations.rotations,
        prim_translations=magnetic_operations.translations,
        prim_time_reversals=magnetic_operations.time_reversals,
    )
    assert magnetic_space_group.uni_number == uni_number
