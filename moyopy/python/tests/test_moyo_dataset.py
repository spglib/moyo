from __future__ import annotations

import numpy as np

import moyopy


def test_moyo_dataset(wurtzite: moyopy.Cell):
    dataset = moyopy.MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480


def test_serialization(wurtzite: moyopy.Cell):
    serialized = wurtzite.serialize_json()
    deserialized = moyopy.Cell.deserialize_json(serialized)
    assert len(wurtzite.positions) == len(deserialized.positions)


def test_moyo_dataset_repr(wurtzite: moyopy.Cell):
    dataset = moyopy.MoyoDataset(wurtzite)
    dataset_str = str(dataset)

    # Test that string representation of MoyoDataset contains key information
    assert "MoyoDataset" in dataset_str
    assert f"number={dataset.number}" in dataset_str
    assert f"hall_number={dataset.hall_number}" in dataset_str
    assert f"operations=<{len(dataset.operations)} operations>" in dataset_str
    assert f"orbits={dataset.orbits}" in dataset_str
    assert f"wyckoffs={dataset.wyckoffs}" in dataset_str

    # Test site_symmetry_symbols content without caring about quote style
    symbols = dataset.site_symmetry_symbols
    assert all(symbol in dataset_str for symbol in symbols)
    assert str(len(symbols)) in dataset_str

    # Test that repr() gives different output
    assert str(dataset) != repr(dataset)


def test_crystal_system(wurtzite: moyopy.Cell):
    # Test wurtzite structure (space group 186, hexagonal)
    # Use higher symprec since wurtzite structure is slightly distorted
    dataset = moyopy.MoyoDataset(wurtzite, symprec=1e-2)
    assert dataset.crystal_system == moyopy.CrystalSystem.Hexagonal
    assert str(dataset.crystal_system) == "Hexagonal"
    assert repr(dataset.crystal_system) == "CrystalSystem.Hexagonal"

    # Test FCC structure (space group 225, cubic)
    positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]
    fcc = moyopy.Cell(basis=np.eye(3), positions=positions, numbers=[0, 0, 0, 0])
    dataset = moyopy.MoyoDataset(fcc)
    crystal_system = dataset.crystal_system
    assert crystal_system == moyopy.CrystalSystem.Cubic

    # Test that crystal_system appears in string representation
    assert f"{crystal_system=!s}" in str(dataset)

    # Test all crystal systems are accessible as enum members
    for crys_sys in (
        "Triclinic",
        "Monoclinic",
        "Orthorhombic",
        "Tetragonal",
        "Trigonal",
        "Hexagonal",
        "Cubic",
    ):
        assert str(getattr(moyopy.CrystalSystem, crys_sys)) == crys_sys

    assert moyopy.CrystalSystem.__module__ == "moyopy"
    assert moyopy.CrystalSystem.__name__ == "CrystalSystem"
    assert (
        repr(moyopy.CrystalSystem) == str(moyopy.CrystalSystem) == "<class 'moyopy.CrystalSystem'>"
    )
