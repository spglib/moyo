from __future__ import annotations

import pytest

from moyopy import PointGroup, operations_from_number


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
