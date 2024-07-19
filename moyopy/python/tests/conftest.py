from math import sqrt

import pytest

import moyopy


@pytest.fixture
def wurtzite() -> moyopy.Cell:
    # https://next-gen.materialsproject.org/materials/mp-560588
    a = 3.81
    c = 6.24
    basis = [
        [a, 0.0, 0.0],
        [-a / 2.0, a * sqrt(3.0) / 2.0, 0.0],
        [0.0, 0.0, c],
    ]
    z1_2b = 0.00014
    z2_2b = 0.37486
    positions = [
        # 2b
        [1 / 3, 2 / 3, z1_2b],
        [2 / 3, 1 / 3, z1_2b + 0.5],
        # 2b
        [1 / 3, 2 / 3, z2_2b],
        [2 / 3, 1 / 3, z2_2b + 0.5],
    ]
    numbers = [1, 1, 2, 2]

    cell = moyopy.Cell(basis, positions, numbers)
    return cell
