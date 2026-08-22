from __future__ import annotations

import pytest

from moyopy import (
    LayerGroup,
    MagneticSpaceGroup,
    PointGroup,
    SpaceGroup,
    UnimodularTransformation,
    integral_normalizer,
    magnetic_operations_from_uni_number,
    operations_from_layer_number,
    operations_from_number,
)


@pytest.mark.parametrize(
    "number,arithmetic_number",
    [
        (158, 45),  # P3c1, 3m1P
        (159, 46),  # P31c, 31mP
    ],
)
def test_identify_point_group(number: int, arithmetic_number: int):
    operations = operations_from_number(number, primitive=True)
    point_group = PointGroup(operations.rotations)
    assert point_group.arithmetic_number == arithmetic_number


@pytest.mark.parametrize(
    "number",
    [
        158,  # P3c1
        159,  # P31c
    ],
)
def test_identify_space_group(number: int):
    operations = operations_from_number(number, primitive=True)
    space_group = SpaceGroup(
        prim_rotations=operations.rotations, prim_translations=operations.translations
    )
    assert space_group.number == number


@pytest.mark.parametrize(
    "number",
    [
        1,  # p1
        8,  # p 2 1 1
        65,  # p 6
    ],
)
def test_identify_layer_group(number: int):
    operations = operations_from_layer_number(number, primitive=True)
    layer_group = LayerGroup(
        prim_rotations=operations.rotations, prim_translations=operations.translations
    )
    assert layer_group.number == number
    assert 1 <= layer_group.hall_number <= 116
    assert len(layer_group.linear) == 3
    assert all(len(row) == 3 for row in layer_group.linear)
    assert len(layer_group.origin_shift) == 3


@pytest.mark.parametrize(
    "number",
    [
        1,  # p1
        8,  # p 2 1 1
    ],
)
def test_identify_layer_group_with_basis(number: int):
    # Pass an explicit identity basis to exercise the from_lattice path.
    # The in-plane Minkowski reduction is the identity for this basis, so
    # the result must agree with the no-basis call.
    operations = operations_from_layer_number(number, primitive=True)
    identity_basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    layer_group = LayerGroup(
        prim_rotations=operations.rotations,
        prim_translations=operations.translations,
        basis=identity_basis,
    )
    assert layer_group.number == number


@pytest.mark.parametrize(
    "uni_number",
    [
        1242,  # R31'_c
    ],
)
def test_identify_magnetic_space_group(uni_number: int):
    magnetic_operations = magnetic_operations_from_uni_number(uni_number, primitive=True)
    magnetic_space_group = MagneticSpaceGroup(
        prim_rotations=magnetic_operations.rotations,
        prim_translations=magnetic_operations.translations,
        prim_time_reversals=magnetic_operations.time_reversals,
    )
    assert magnetic_space_group.uni_number == uni_number


@pytest.mark.parametrize(
    "number,expected",
    [
        pytest.param(16, 6, id="P 2 2 2"),  # Permutations for the three rotation axes
        pytest.param(17, 2, id="P 2 2 2_1"),  # Permutations for the two rotation axes
        pytest.param(18, 2, id="P 2_1 2_1 2"),  # Permutations for the two skew rotation axes
        pytest.param(19, 6, id="P 2_1 2_1 2_1"),  # Permutations for the three skew rotation axes
    ],
)
def test_integral_normalizer_orthorhombic(number: int, expected: int):
    operations = operations_from_number(number, primitive=True)
    normalizer = integral_normalizer(operations.rotations, operations.translations)
    assert len(normalizer) == expected


def test_integral_normalizer_defaults_to_small_generators():
    operations = operations_from_number(158, primitive=True)

    actual = integral_normalizer(
        operations.rotations,
        operations.translations,
    )
    expected = integral_normalizer(
        operations.rotations,
        operations.translations,
        prim_generators=list(range(len(operations))),
    )

    assert actual
    assert {transformation.serialize_json() for transformation in actual} == {
        transformation.serialize_json() for transformation in expected
    }
    assert all(isinstance(transformation, UnimodularTransformation) for transformation in actual)
    assert all(len(transformation.linear) == 3 for transformation in actual)
    assert all(len(row) == 3 for transformation in actual for row in transformation.linear)
    assert all(len(transformation.origin_shift) == 3 for transformation in actual)


def test_integral_normalizer_invalid_generator_index():
    operations = operations_from_number(158, primitive=True)

    with pytest.raises(ValueError, match="out of range"):
        integral_normalizer(
            operations.rotations,
            operations.translations,
            prim_generators=[len(operations)],
        )
