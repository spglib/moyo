from __future__ import annotations

import numpy as np
import pytest

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


def test_from_py_obj(wurtzite: moyopy.Cell):
    # Test with pymatgen Structure
    structure = MoyoAdapter.get_structure(wurtzite)
    cell1 = MoyoAdapter.from_py_obj(structure)
    assert np.allclose(cell1.basis, wurtzite.basis)
    assert np.allclose(cell1.positions, wurtzite.positions)
    assert cell1.numbers == wurtzite.numbers

    # Test with ASE Atoms
    atoms = MoyoAdapter.get_atoms(wurtzite)
    cell2 = MoyoAdapter.from_py_obj(atoms)
    assert np.allclose(cell2.basis, wurtzite.basis)
    assert np.allclose(cell2.positions, wurtzite.positions)
    assert cell2.numbers == wurtzite.numbers

    # Test with MSONAtoms
    mson_atoms = structure.to_ase_atoms()
    cell3 = MoyoAdapter.from_py_obj(mson_atoms)
    assert np.allclose(cell3.basis, wurtzite.basis)
    assert np.allclose(cell3.positions, wurtzite.positions)
    assert cell3.numbers == wurtzite.numbers

    # Test invalid input type
    with pytest.raises(TypeError, match="Expected Structure or Atoms"):
        MoyoAdapter.from_py_obj([1, 2, 3])
