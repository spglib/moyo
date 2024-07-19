from __future__ import annotations

import numpy as np

import moyopy
from moyopy.interface import MoyoAdapter


def test_pymatgen_moyopy(wurtzite: moyopy.Cell):
    structure = MoyoAdapter.get_structure(wurtzite)
    cell = MoyoAdapter.from_structure(structure)
    assert np.allclose(cell.basis, wurtzite.basis)
    assert np.allclose(cell.positions, wurtzite.positions)
    assert cell.numbers == wurtzite.numbers


def test_ase_moyopy(wurtzite: moyopy.Cell):
    atoms = MoyoAdapter.get_atoms(wurtzite)
    cell = MoyoAdapter.from_atoms(atoms)
    assert np.allclose(cell.basis, wurtzite.basis)
    assert np.allclose(cell.positions, wurtzite.positions)
    assert cell.numbers == wurtzite.numbers
