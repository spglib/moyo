from __future__ import annotations

from moyopy import operations_from_number


def test_operations_from_number():
    operations = operations_from_number(number=230)  # Ia-3d
    assert operations.num_operations == 48 * 2
