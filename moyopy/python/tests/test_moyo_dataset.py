from __future__ import annotations

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
    repr_str = repr(dataset)

    # Test that repr contains key information
    assert "MoyoDataset" in repr_str
    assert f"number={dataset.number}" in repr_str
    assert f"hall_number={dataset.hall_number}" in repr_str
    assert f"operations=<{len(dataset.operations)} operations>" in repr_str
    assert f"orbits={dataset.orbits}" in repr_str
    assert f"wyckoffs={dataset.wyckoffs}" in repr_str

    # Test site_symmetry_symbols content without caring about quote style
    symbols = dataset.site_symmetry_symbols
    assert all(symbol in repr_str for symbol in symbols)
    assert str(len(symbols)) in repr_str

    # Test that str() gives same output as repr()
    assert str(dataset) == repr_str
