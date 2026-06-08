from __future__ import annotations

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
