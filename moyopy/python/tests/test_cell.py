from moyopy import Cell


def test_cell_serialization(wurtzite: Cell):
    serialized = wurtzite.serialize_json()
    deserialized = Cell.deserialize_json(serialized)
    assert len(wurtzite.positions) == len(deserialized.positions)


def test_cell_py_obj_serialization(wurtzite: Cell):
    deserialized = wurtzite.as_dict()
    serialized = Cell.from_dict(deserialized)
    assert len(wurtzite.positions) == len(serialized.positions)
