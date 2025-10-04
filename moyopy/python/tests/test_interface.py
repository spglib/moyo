from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pymatgen.core import Composition, DummySpecies, Element, Species, Structure

from moyopy.interface import MoyoAdapter, MoyoNonCollinearMagneticAdapter

if TYPE_CHECKING:
    from moyopy import Cell, NonCollinearMagneticCell


def test_pymatgen_moyopy(wurtzite: Cell):
    structure = MoyoAdapter.get_structure(wurtzite)
    cell = MoyoAdapter.from_structure(structure)
    assert np.allclose(cell.basis, wurtzite.basis)
    assert np.allclose(cell.positions, wurtzite.positions)
    assert cell.numbers == wurtzite.numbers


def test_pymatgen_disordered_moyopy():
    structure = Structure(
        lattice=4.0 * np.eye(3),
        coords=[
            [0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ],
        species=[
            Element("Cu"),
            Species("Cu2+"),
            DummySpecies("X0+"),
            Composition("Cu0.5Au0.5"),
        ],
    )
    cell, unique_species_mapping = MoyoAdapter.from_disordered_structure(structure)
    structure2 = MoyoAdapter.get_structure(cell, unique_species_mapping=unique_species_mapping)
    assert np.allclose(structure.lattice.matrix, structure2.lattice.matrix)
    assert np.allclose(structure.frac_coords, structure2.frac_coords)
    assert structure.species_and_occu == structure2.species_and_occu


def test_ase_moyopy(wurtzite: Cell):
    atoms = MoyoAdapter.get_atoms(wurtzite)
    cell = MoyoAdapter.from_atoms(atoms)
    assert np.allclose(cell.basis, wurtzite.basis)
    assert np.allclose(cell.positions, wurtzite.positions)
    assert cell.numbers == wurtzite.numbers


def test_from_py_obj(wurtzite: Cell):
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


def test_noncollinear_adapter(rutile_type3_noncollinear: NonCollinearMagneticCell):
    structure = MoyoNonCollinearMagneticAdapter.get_structure(rutile_type3_noncollinear)
    magnetic_cell = MoyoNonCollinearMagneticAdapter.from_structure(structure)
    assert np.allclose(magnetic_cell.basis, rutile_type3_noncollinear.basis)
    assert np.allclose(magnetic_cell.positions, rutile_type3_noncollinear.positions)
    assert magnetic_cell.numbers == rutile_type3_noncollinear.numbers
    assert np.allclose(magnetic_cell.magnetic_moments, rutile_type3_noncollinear.magnetic_moments)
