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
