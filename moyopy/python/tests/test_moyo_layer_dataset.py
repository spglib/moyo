from __future__ import annotations

from math import cos, radians, sin

import pytest

from moyopy import Cell, LayerSetting, MoyoLayerDataset


@pytest.fixture
def p1_layer_cell() -> Cell:
    # Two unrelated atoms in an oblique 2D lattice with c perpendicular.
    a = 1.0
    b = 1.5
    gamma = radians(80.0)
    basis = [
        [a, 0.0, 0.0],
        [b * cos(gamma), b * sin(gamma), 0.0],
        [0.0, 0.0, 5.0],
    ]
    positions = [
        [0.1, 0.2, 0.1],
        [0.3, 0.5, 0.2],
    ]
    numbers = [1, 1]
    return Cell(basis, positions, numbers)


@pytest.fixture
def p4mm_layer_cell() -> Cell:
    # Single atom at the 1a site of LG 55 (p4mm).
    basis = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 5.0],
    ]
    positions = [[0.0, 0.0, 0.1]]
    numbers = [1]
    return Cell(basis, positions, numbers)


def test_moyo_layer_dataset_p1(p1_layer_cell: Cell):
    # The two atoms accidentally land at inversion-related general positions
    # through (0.2, 0.35, 0.15), so the layer group is LG 2 (p-1) rather
    # than LG 1; the off-origin inversion has a non-zero free t_z.
    dataset = MoyoLayerDataset(p1_layer_cell)
    assert dataset.number == 2
    assert dataset.hall_number == 2
    assert dataset.pearson_symbol == "mp2"
    assert dataset.orbits == [0, 0]
    assert dataset.wyckoffs == ["e", "e"]
    # |c| preserved and c along z in std_cell.
    std_basis = dataset.std_cell.basis
    assert std_basis[2][2] == pytest.approx(5.0)
    assert std_basis[2][0] == pytest.approx(0.0)
    assert std_basis[2][1] == pytest.approx(0.0)


def test_moyo_layer_dataset_p4mm(p4mm_layer_cell: Cell):
    # A single atom at (0, 0, 0.1) sits on a horizontal mirror plane, so
    # the layer group is the centro-symmetric supergroup LG 61 (p4/mmm),
    # not LG 55 (p4mm).
    dataset = MoyoLayerDataset(p4mm_layer_cell, setting=LayerSetting.standard())
    assert dataset.number == 61
    assert dataset.pearson_symbol.startswith("tp")


def test_moyo_layer_dataset_serialization(p1_layer_cell: Cell):
    dataset = MoyoLayerDataset(p1_layer_cell)
    serialized = dataset.serialize_json()
    deserialized = MoyoLayerDataset.deserialize_json(serialized)
    assert deserialized.number == dataset.number
    assert deserialized.hall_number == dataset.hall_number
    assert deserialized.std_cell.num_atoms == dataset.std_cell.num_atoms


def test_moyo_layer_dataset_as_dict_round_trip(p1_layer_cell: Cell):
    dataset = MoyoLayerDataset(p1_layer_cell)
    obj = dataset.as_dict()
    restored = MoyoLayerDataset.from_dict(obj)
    assert restored.number == dataset.number
    assert restored.hall_number == dataset.hall_number


def test_moyo_layer_dataset_repr(p1_layer_cell: Cell):
    dataset = MoyoLayerDataset(p1_layer_cell)
    s = str(dataset)
    assert "MoyoLayerDataset" in s
    assert f"number={dataset.number}" in s
    assert f"hall_number={dataset.hall_number}" in s


def test_layer_setting_variants():
    spglib = LayerSetting.spglib()
    standard = LayerSetting.standard()
    explicit = LayerSetting.hall_number(1)
    for s in (spglib, standard, explicit):
        assert isinstance(s.serialize_json(), str)


def test_moyo_layer_dataset_rejects_tilted_c():
    # c with a non-zero in-plane component must be rejected.
    basis = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.0, 5.0],
    ]
    cell = Cell(basis, [[0.0, 0.0, 0.0]], [1])
    with pytest.raises(ValueError):
        MoyoLayerDataset(cell)
