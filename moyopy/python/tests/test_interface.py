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
    # Create different structure types
    structure = MoyoAdapter.get_structure(wurtzite)
    atoms = MoyoAdapter.get_atoms(wurtzite)
    mson_atoms = structure.to_ase_atoms()

    # Test conversion from each type
    for struct in [structure, atoms, mson_atoms]:
        cell = MoyoAdapter.from_py_obj(struct)
        assert np.allclose(cell.basis, wurtzite.basis)
        assert np.allclose(cell.positions, wurtzite.positions)
        assert cell.numbers == wurtzite.numbers

    # Test invalid input type
    with pytest.raises(TypeError, match="Expected Structure, Atoms, or MSONAtoms, got list"):
        MoyoAdapter.from_py_obj([1, 2, 3])
