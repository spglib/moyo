from __future__ import annotations

import logging

import numpy as np
import pytest

from moyopy import (
    Cell,
    CollinearMagneticCell,
    MoyoCollinearMagneticDataset,
    MoyoDataset,
    NormalizerWyckoffPositions,
    UnimodularTransformation,
    WyckoffPosition,
)


def _word(wyckoffs: list[WyckoffPosition]) -> str:
    return "".join(w.letter for w in wyckoffs)


def test_moyo_dataset(wurtzite: Cell):
    dataset = MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480
    assert dataset.pearson_symbol == "hP4"


@pytest.mark.parametrize("n", [1, 3], ids=["conventional", "3x1x1"])
@pytest.mark.parametrize("magnetic", [False, True], ids=["nonmagnetic", "magnetic"])
def test_supercell_rotation_warning(caplog: pytest.LogCaptureFixture, n: int, magnetic: bool):
    base = [
        ([0.0, 0.0, 0.0], 11),
        ([0.5, 0.5, 0.0], 11),
        ([0.5, 0.0, 0.5], 11),
        ([0.0, 0.5, 0.5], 11),
        ([0.5, 0.5, 0.5], 17),
        ([0.0, 0.0, 0.5], 17),
        ([0.0, 0.5, 0.0], 17),
        ([0.5, 0.0, 0.0], 17),
    ]
    basis = [[5.64 * n, 0.0, 0.0], [0.0, 5.64, 0.0], [0.0, 0.0, 5.64]]
    positions = [[(p[0] + copy) / n, p[1], p[2]] for copy in range(n) for p, _ in base]
    numbers = [number for _ in range(n) for _, number in base]
    logger = "moyo"

    with caplog.at_level(logging.WARNING, logger=logger):
        if magnetic:
            cell = CollinearMagneticCell(basis, positions, numbers, [0.0] * len(numbers))
            operations = MoyoCollinearMagneticDataset(cell, symprec=1e-3).magnetic_operations
            # Zero moments allow both time-reversal states for every spatial operation.
            assert len(operations) == 384
        else:
            operations = MoyoDataset(Cell(basis, positions, numbers), symprec=1e-3).operations
            assert len(operations) == 192

    records = [record for record in caplog.records if record.name == logger]
    if n == 1:
        assert records == []
    else:
        (record,) = records
        assert record.levelno == logging.WARNING
        assert "non-integer rotation matrices in the input-cell basis" in record.getMessage()
        assert (
            "returning only operations compatible with the input-cell lattice"
            in record.getMessage()
        )
        if not magnetic:
            assert "orbits use primitive-cell symmetry" in record.getMessage()
            assert "may differ from spglib's equivalent_atoms" in record.getMessage()

    metric = np.array(basis) @ np.array(basis).T
    for rotation in np.array(operations.rotations):
        np.testing.assert_allclose(rotation.T @ metric @ rotation, metric, atol=1e-8)


@pytest.mark.parametrize("ny, expected_operations", [(1, 32), (3, 48)])
def test_antiferromagnetic_supercell_warning(
    caplog: pytest.LogCaptureFixture, ny: int, expected_operations: int
):
    # Alternating moments double the primitive magnetic cell along x. Preparing
    # candidates from the nonmagnetic cell must not emit an intermediate warning.
    cell = CollinearMagneticCell(
        [[2.0, 0.0, 0.0], [0.0, float(ny), 0.0], [0.0, 0.0, 1.0]],
        [[i / 2, j / ny, 0.0] for i in range(2) for j in range(ny)],
        [1] * (2 * ny),
        [(-1.0) ** i for i in range(2) for _ in range(ny)],
    )
    with caplog.at_level(logging.WARNING, logger="moyo"):
        dataset = MoyoCollinearMagneticDataset(cell)

    assert len(dataset.magnetic_operations) == expected_operations
    records = [record for record in caplog.records if record.name.startswith("moyo")]
    assert len(records) == (ny > 1)
    for record in records:
        assert record.levelno == logging.WARNING
        assert "magnetic symmetry operations with non-integer rotation" in record.getMessage()


def test_moyo_dataset_serialization(wurtzite: Cell):
    dataset = MoyoDataset(wurtzite)
    serialized = dataset.serialize_json()
    deserialized = MoyoDataset.deserialize_json(serialized)
    assert deserialized.number == dataset.number
    assert deserialized.std_cell.num_atoms == dataset.std_cell.num_atoms


def test_moyo_dataset_py_obj_serialization(wurtzite: Cell):
    dataset = MoyoDataset(wurtzite)
    deserialized = dataset.as_dict()
    serialized = MoyoDataset.from_dict(deserialized)
    assert serialized.number == dataset.number
    assert serialized.std_cell.num_atoms == dataset.std_cell.num_atoms


def test_moyo_dataset_repr(wurtzite: Cell):
    dataset = MoyoDataset(wurtzite)
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


def test_normalizer_wyckoff_positions_perovskite():
    # Cubic perovskite ABO3 (Pm-3m, #221): A at 1a, B at 1b, O at 3c.
    a = 4.0
    basis = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]
    positions = [
        [0.0, 0.0, 0.0],  # A (1a)
        [0.5, 0.5, 0.5],  # B (1b)
        [0.5, 0.5, 0.0],  # O (3c)
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ]
    numbers = [0, 1, 2, 2, 2]
    dataset = MoyoDataset(Cell(basis, positions, numbers))
    assert dataset.number == 221
    assert dataset.wyckoffs == ["a", "b", "c", "c", "c"]

    result = dataset.normalizer_wyckoff_positions()
    assert isinstance(result, NormalizerWyckoffPositions)

    # Identity setting reproduces the dataset's own Wyckoff sequence.
    assert _word(result.wyckoffs) == "abccc"
    assert [w.multiplicity for w in result.wyckoffs] == [1, 1, 3, 3, 3]

    # The (1/2,1/2,1/2) normalizer translation swaps a<->b and c<->d, so there
    # are exactly two distinct sequences.
    sequences = {_word(seq) for _, seq in result.coset_representatives}
    assert sequences == {"abccc", "baddd"}
    assert len(result.coset_representatives) == 2

    # First coset representative is the identity paired with `wyckoffs`.
    op0, seq0 = result.coset_representatives[0]
    assert isinstance(op0, UnimodularTransformation)
    assert _word(seq0) == "abccc"

    # The sequences and stabilizer size are basis-independent.
    result_prim = dataset.normalizer_wyckoff_positions(primitive=True)
    assert {_word(seq) for _, seq in result_prim.coset_representatives} == sequences
    assert len(result_prim.stabilizer) == len(result.stabilizer)


def test_moyo_collinear_magnetic_dataset(rutile_type3: CollinearMagneticCell):
    uni_number = 1158
    dataset = MoyoCollinearMagneticDataset(rutile_type3)
    assert dataset.uni_number == uni_number
    dataset2 = MoyoCollinearMagneticDataset(dataset.std_mag_cell)
    assert dataset2.uni_number == uni_number
    dataset3 = MoyoCollinearMagneticDataset(dataset.prim_std_mag_cell)
    assert dataset3.uni_number == uni_number
