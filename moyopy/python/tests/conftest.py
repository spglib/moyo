from math import sqrt

import pytest

from moyopy import Cell, CollinearMagneticCell


@pytest.fixture
def wurtzite() -> Cell:
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

    cell = Cell(basis, positions, numbers)
    return cell


@pytest.fixture
def rutile_type3() -> CollinearMagneticCell:
    basis = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    positions = [
        # Ti (2a)
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        # O (4f)
        [0.3, 0.3, 0.0],
        [0.7, 0.7, 0.0],
        [0.2, 0.8, 0.5],
        [0.8, 0.2, 0.5],
    ]
    numbers = [0, 0, 1, 1, 1, 1]
    magnetic_moments = [0.7, -0.7, 0.0, 0.0, 0.0, 0.0]

    magnetic_cell = CollinearMagneticCell(basis, positions, numbers, magnetic_moments)
    return magnetic_cell
