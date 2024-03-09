from __future__ import annotations

from math import sqrt

import pytest

import moyo


@pytest.fixture
def wurtzite() -> moyo.Structure:
    # https://next-gen.materialsproject.org/materials/mp-560588
    a = 3.81
    c = 6.24
    basis = [
        [a, 0.0, 0.0],
        [-a / 2.0, a * sqrt(3.0) / 2.0, 0.0],
        [0.0, 0.0, c],
    ]
    positions = [
        [1 / 3, 2 / 3, 0.00014],
        [1 / 3, 2 / 3, 0.37486],
    ]
    numbers = [0, 1]

    structure = moyo.Structure(basis, positions, numbers)
    return structure
