from __future__ import annotations

from moyopy import Cell, CollinearMagneticCell, MoyoCollinearMagneticDataset, MoyoDataset


def test_moyo_dataset(wurtzite: Cell):
    dataset = MoyoDataset(wurtzite)
    assert dataset.number == 186
    assert dataset.hall_number == 480
    assert dataset.pearson_symbol == "hP4"


def test_serialization(wurtzite: Cell):
    serialized = wurtzite.serialize_json()
    deserialized = Cell.deserialize_json(serialized)
    assert len(wurtzite.positions) == len(deserialized.positions)


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


def test_moyo_collinear_magnetic_dataset(rutile_type3: CollinearMagneticCell):
    dataset = MoyoCollinearMagneticDataset(rutile_type3)
    assert dataset.uni_number == 1158
