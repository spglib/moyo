from __future__ import annotations

import numpy as np
import pytest
from _helpers import REBASES, assert_proper_rotation

from moyopy import Cell, MoyoDataset, Setting, SpaceGroup, operations_from_number

ENANTIOMORPHIC_PAIRS = [
    (76, 78),
    (91, 95),
    (92, 96),
    (144, 145),
    (151, 153),
    (152, 154),
    (169, 170),
    (171, 172),
    (178, 179),
    (180, 181),
    (212, 213),
]
# Fixed database settings for the fixture groups, independent of identification.
HALL_NUMBERS = {
    5: 9,
    76: 350,
    78: 352,
    91: 368,
    92: 369,
    95: 372,
    96: 373,
    144: 431,
    145: 432,
    151: 440,
    152: 441,
    153: 442,
    154: 443,
    169: 463,
    170: 464,
    171: 465,
    172: 466,
    178: 472,
    179: 473,
    180: 474,
    181: 475,
    212: 508,
    213: 509,
}


def _operation_keys(rotations, translations):
    return {
        (tuple(np.asarray(w).flatten()), tuple(np.rint(np.asarray(t) * 1e8).astype(int) % 10**8))
        for w, t in zip(rotations, translations)
    }


def _transform_operations(rotations, translations, linear, origin):
    inverse = np.linalg.inv(linear)
    rotations = np.asarray(rotations)
    translations = np.asarray(translations)
    transformed_rotations = inverse @ rotations @ linear
    np.testing.assert_allclose(transformed_rotations, np.rint(transformed_rotations), atol=1e-8)
    transformed_translations = (rotations @ origin + translations - origin) @ inverse.T
    return np.rint(transformed_rotations).astype(int), transformed_translations


@pytest.mark.parametrize("pair", ENANTIOMORPHIC_PAIRS, ids=lambda pair: f"{pair[0]}-{pair[1]}")
@pytest.mark.parametrize("mirrored", [False, True], ids=["original", "mirror"])
@pytest.mark.parametrize("rebase", REBASES)
def test_enantiomorphic_handedness(pair, mirrored, rebase):
    number, mirror_number = pair
    expected_number = mirror_number if mirrored else number
    if number >= 195:
        basis = np.diag([4.0, 4.0, 4.0])
    elif number >= 143:
        basis = np.array([[4.0, -2.0, 0.0], [0.0, np.sqrt(12), 0.0], [0.0, 0.0, 6.0]])
    else:
        basis = np.diag([4.0, 4.0, 6.0])

    reference = operations_from_number(number, primitive=True)
    rotations = np.array(reference.rotations)
    translations = np.array(reference.translations)
    # Two generic, differently labelled orbits avoid accidental extra symmetries.
    positions = np.concatenate(
        [
            rotations @ site + translations
            for site in ([0.137, 0.271, 0.389], [0.219, 0.413, 0.157])
        ]
    )
    numbers = np.repeat([1, 2], len(reference))
    # This is an active Cartesian reflection, so it must exchange the pair.
    if mirrored:
        basis = np.diag([-1, 1, 1]) @ basis
    cartesian_rotation = np.array([[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]])
    basis = cartesian_rotation @ basis

    origin = np.array([0.13, 0.27, 0.19])
    input_basis = basis @ rebase
    input_positions = (positions - origin) @ np.linalg.inv(rebase).T
    # Rebasing and shifting the origin preserve the Cartesian sites.
    np.testing.assert_allclose(
        input_positions @ input_basis.T + basis @ origin, positions @ basis.T, atol=1e-12
    )
    input_positions %= 1
    cell = Cell(input_basis.T.tolist(), input_positions.tolist(), numbers.tolist())
    expected_hall_number = HALL_NUMBERS[expected_number]
    input_rotations, input_translations = _transform_operations(
        rotations, translations, rebase, origin
    )
    target = operations_from_number(expected_number, primitive=True)
    target_keys = _operation_keys(target.rotations, target.translations)
    sign = np.sign(np.linalg.det(input_basis))

    for setting in [None, Setting.spglib(), Setting.hall_number(expected_hall_number)]:
        group = SpaceGroup(
            input_rotations.tolist(),
            input_translations.tolist(),
            basis=cell.basis,
            setting=setting,
        )
        assert group.number == expected_number
        assert group.hall_number == expected_hall_number
        assert round(np.linalg.det(group.linear)) == sign
        assert (
            _operation_keys(
                *_transform_operations(
                    input_rotations, input_translations, group.linear, np.array(group.origin_shift)
                )
            )
            == target_keys
        )

        # Identification is independent of rotate_basis; only datasets need both modes.
        for rotate_basis in [False, True]:
            try:
                dataset = MoyoDataset(
                    cell, setting=setting, rotate_basis=rotate_basis, symprec=1e-5
                )
                assert dataset.number == expected_number
                assert dataset.hall_number == expected_hall_number
                assert _operation_keys(
                    dataset.operations.rotations, dataset.operations.translations
                ) == (_operation_keys(input_rotations, input_translations))
                # Compare partitions, rather than arbitrary representative site indices.
                same_orbit = np.equal.outer(dataset.orbits, dataset.orbits)
                np.testing.assert_array_equal(same_orbit, np.equal.outer(numbers, numbers))
                rotation = np.array(dataset.std_rotation_matrix)
                assert_proper_rotation(rotation, rotate_basis)

                for standardized, linear, shift in [
                    (dataset.std_cell, dataset.std_linear, dataset.std_origin_shift),
                    (
                        dataset.prim_std_cell,
                        dataset.prim_std_linear,
                        dataset.prim_std_origin_shift,
                    ),
                ]:
                    assert np.linalg.det(standardized.basis) > 0
                    assert np.sign(np.linalg.det(linear)) == sign
                    np.testing.assert_allclose(
                        np.array(standardized.basis).T, rotation @ input_basis @ linear, atol=1e-8
                    )
                    transformed_positions = (input_positions - shift) @ np.linalg.inv(linear).T
                    differences = transformed_positions[:, None] - np.array(standardized.positions)
                    differences -= np.rint(differences)
                    matches = (np.linalg.norm(differences, axis=-1) < 1e-7) & (
                        numbers[:, None] == np.array(standardized.numbers)
                    )
                    assert np.all(matches.sum(axis=1) == 1)
                    assert (
                        _operation_keys(
                            *_transform_operations(
                                input_rotations, input_translations, linear, np.array(shift)
                            )
                        )
                        == target_keys
                    )

                mapping = np.array(dataset.mapping_std_prim)
                prim_positions = np.array(dataset.prim_std_cell.positions)[mapping]
                expected_positions = (
                    input_positions - dataset.prim_std_origin_shift
                ) @ np.linalg.inv(dataset.prim_std_linear).T
                diff = prim_positions - expected_positions
                np.testing.assert_allclose(diff - np.rint(diff), 0, atol=1e-7)
                np.testing.assert_array_equal(
                    np.array(dataset.prim_std_cell.numbers)[mapping], numbers
                )
            except Exception as error:
                raise AssertionError(f"setting={setting}, rotate_basis={rotate_basis}") from error


@pytest.mark.parametrize("number,index,centering", [(5, 1, 2), (76, 2, 1)])
@pytest.mark.parametrize("rotate_basis", [False, True])
def test_centered_and_supercell_handedness(number, index, centering, rotate_basis):
    reference = operations_from_number(number)
    rotations = np.array(reference.rotations)
    translations = np.array(reference.translations)
    positions = np.concatenate(
        [
            rotations @ site + translations
            for site in ([0.137, 0.271, 0.389], [0.219, 0.413, 0.157])
        ]
    )
    numbers = np.repeat([1, 2], len(reference))
    basis = (
        np.array([[3.0, 0.0, -1.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]])
        if number == 5
        else np.diag([4.0, 4.0, 6.0])
    )
    supercell = np.diag([1, 1, index])
    positions = np.concatenate(
        [(positions + [0, 0, k]) @ np.linalg.inv(supercell).T for k in range(index)]
    )
    numbers = np.tile(numbers, index)
    # Translation equivalence is known from the original primitive lattice.
    # For C centering, (a+b)/2, (-a+b)/2, c form a primitive basis.
    to_primitive = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]) if centering == 2 else np.eye(3)
    primitive_positions = positions @ supercell.T @ to_primitive.T
    differences = primitive_positions[:, None] - primitive_positions
    differences -= np.rint(differences)
    same_primitive_site = (np.linalg.norm(differences, axis=-1) < 1e-8) & np.equal.outer(
        numbers, numbers
    )
    assert np.all(same_primitive_site.sum(axis=1) == index * centering)
    # Lift the known operations through the supercell's translation cosets.
    supercell_rotations, supercell_translations = _transform_operations(
        np.tile(rotations, (index, 1, 1)),
        np.concatenate([translations + [0, 0, k] for k in range(index)]),
        supercell,
        np.zeros(3),
    )
    basis = basis @ supercell
    rebase = np.array([[0, -1, 2], [1, 0, 1], [0, 0, -1]])
    origin = np.array([0.13, 0.27, 0.19])
    positions = (positions - origin) @ np.linalg.inv(rebase).T
    input_rotations, input_translations = _transform_operations(
        supercell_rotations, supercell_translations, rebase, origin
    )
    input_basis = basis @ rebase
    cell = Cell(input_basis.T.tolist(), positions.tolist(), numbers.tolist())
    dataset = MoyoDataset(cell, symprec=1e-5, rotate_basis=rotate_basis)
    assert dataset.number == number
    assert dataset.hall_number == HALL_NUMBERS[number]
    assert _operation_keys(dataset.operations.rotations, dataset.operations.translations) == (
        _operation_keys(input_rotations, input_translations)
    )
    np.testing.assert_array_equal(
        np.equal.outer(dataset.orbits, dataset.orbits), np.equal.outer(numbers, numbers)
    )
    np.testing.assert_array_equal(
        np.equal.outer(dataset.mapping_std_prim, dataset.mapping_std_prim),
        same_primitive_site,
    )
    for standardized, linear, expected_volume in [
        (dataset.std_cell, dataset.std_linear, abs(np.linalg.det(basis)) / index),
        (
            dataset.prim_std_cell,
            dataset.prim_std_linear,
            abs(np.linalg.det(basis)) / index / centering,
        ),
    ]:
        assert np.linalg.det(standardized.basis) == pytest.approx(expected_volume)
        assert np.linalg.det(linear) < 0
        np.testing.assert_allclose(
            np.array(standardized.basis).T,
            np.array(dataset.std_rotation_matrix) @ input_basis @ linear,
            atol=1e-8,
        )
    assert np.linalg.det(dataset.std_rotation_matrix) == pytest.approx(1)
    prim_positions = np.array(dataset.prim_std_cell.positions)[dataset.mapping_std_prim]
    expected_positions = (positions - dataset.prim_std_origin_shift) @ np.linalg.inv(
        dataset.prim_std_linear
    ).T
    diff = prim_positions - expected_positions
    np.testing.assert_allclose(diff - np.rint(diff), 0, atol=1e-7)
    np.testing.assert_array_equal(
        np.array(dataset.prim_std_cell.numbers)[dataset.mapping_std_prim], numbers
    )


def test_operations_only_retains_reference_orientation():
    # Without a lattice, the same fractional operations always use a right-handed
    # reference orientation. They cannot reveal the actual Cartesian handedness.
    reference = operations_from_number(76, primitive=True)
    group = SpaceGroup(reference.rotations, reference.translations)
    assert group.number == 76
    assert round(np.linalg.det(group.linear)) == 1
    rotations, translations = _transform_operations(
        reference.rotations, reference.translations, np.diag([-1, 1, 1]), np.zeros(3)
    )
    reflected = SpaceGroup(rotations.tolist(), translations.tolist())
    assert reflected.number == 78
    assert round(np.linalg.det(reflected.linear)) == 1


@pytest.mark.parametrize("number,hall_number,centering", [(136, 419, 1), (225, 523, 4)])
@pytest.mark.parametrize("perturbed", [False, True], ids=["exact", "perturbed"])
@pytest.mark.parametrize("rebase", REBASES)
@pytest.mark.parametrize("rotate_basis", [False, True])
def test_refined_wyckoff_handedness(
    number, hall_number, centering, perturbed, rebase, rotate_basis
):
    if number == 136:
        # Rutile: 2a + 4f in P4_2/mnm, including a free Wyckoff parameter.
        u = 0.3
        basis = np.diag([4.6, 4.6, 2.95])
        positions = np.array(
            [
                [0, 0, 0],
                [0.5, 0.5, 0.5],
                [u, u, 0],
                [-u, -u, 0],
                [0.5 + u, 0.5 - u, 0.5],
                [0.5 - u, 0.5 + u, 0.5],
            ]
        )
        numbers = np.array([1, 1, 2, 2, 2, 2])
        site_symmetries = {1: "m.mm", 2: "m.2m"}
    else:
        # Rocksalt: two special orbits and an explicit F-centered primitive map.
        basis = np.eye(3) * 5.6
        face_centers = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
        positions = np.concatenate([face_centers, face_centers + 0.5])
        numbers = np.repeat([1, 2], 4)
        site_symmetries = {1: "m-3m", 2: "m-3m"}

    to_primitive = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) if centering == 4 else np.eye(3)
    primitive_positions = positions @ to_primitive.T
    differences = primitive_positions[:, None] - primitive_positions
    differences -= np.rint(differences)
    same_primitive_site = (np.linalg.norm(differences, axis=-1) < 1e-8) & np.equal.outer(
        numbers, numbers
    )
    if perturbed:
        basis = (
            np.array(
                [[1.0007, 0.0002, -0.0003], [0.0004, 0.9995, 0.0001], [0.0002, -0.0003, 1.0003]]
            )
            @ basis
        )
        positions += 1e-5 * np.sin(np.arange(positions.size).reshape(positions.shape))
    cartesian_rotation = np.array([[1, 0, 0], [0, 0.6, -0.8], [0, 0.8, 0.6]])
    input_basis = cartesian_rotation @ basis @ rebase
    input_positions = (positions - [0.13, 0.27, 0.19]) @ np.linalg.inv(rebase).T
    dataset = MoyoDataset(
        Cell(input_basis.T.tolist(), input_positions.tolist(), numbers.tolist()),
        symprec=1e-2,
        rotate_basis=rotate_basis,
    )
    assert (dataset.number, dataset.hall_number) == (number, hall_number)
    np.testing.assert_array_equal(
        np.equal.outer(dataset.orbits, dataset.orbits), np.equal.outer(numbers, numbers)
    )
    np.testing.assert_array_equal(
        np.equal.outer(dataset.mapping_std_prim, dataset.mapping_std_prim), same_primitive_site
    )
    rotation = np.array(dataset.std_rotation_matrix)
    assert_proper_rotation(rotation, rotate_basis)
    for cell, linear in [
        (dataset.std_cell, dataset.std_linear),
        (dataset.prim_std_cell, dataset.prim_std_linear),
    ]:
        assert np.linalg.det(cell.basis) > 0
        assert np.sign(np.linalg.det(linear)) == np.sign(np.linalg.det(input_basis))
        stretch = rotation.T @ np.array(cell.basis).T @ np.linalg.inv(input_basis @ linear)
        np.testing.assert_allclose(stretch, stretch.T, atol=1e-8)
        assert np.all(np.linalg.eigvalsh(stretch) > 0)
        if perturbed:
            assert np.linalg.norm(stretch - np.eye(3)) > 1e-4
        else:
            np.testing.assert_allclose(stretch, np.eye(3), atol=1e-8)
    np.testing.assert_allclose(
        np.array(dataset.prim_std_cell.basis).T @ to_primitive,
        np.array(dataset.std_cell.basis).T,
        atol=1e-8,
    )
    mapped = np.array(dataset.prim_std_cell.positions)[dataset.mapping_std_prim]
    selected = (input_positions - dataset.prim_std_origin_shift) @ np.linalg.inv(
        dataset.prim_std_linear
    ).T
    residuals = mapped - selected
    residuals -= np.rint(residuals)
    residual_norms = np.linalg.norm(residuals @ np.array(dataset.prim_std_cell.basis), axis=1)
    assert np.max(residual_norms) < (1e-3 if perturbed else 1e-8)
    if perturbed:
        assert np.max(residual_norms) > 1e-6
    np.testing.assert_array_equal(
        np.array(dataset.prim_std_cell.numbers)[dataset.mapping_std_prim], numbers
    )

    reference = operations_from_number(number)
    rotations, translations = np.array(reference.rotations), np.array(reference.translations)
    std_basis = np.array(dataset.std_cell.basis).T
    metric = std_basis.T @ std_basis
    std_positions, std_numbers = (
        np.array(dataset.std_cell.positions),
        np.array(dataset.std_cell.numbers),
    )
    for w, t in zip(rotations, translations):
        np.testing.assert_allclose(w.T @ metric @ w, metric, atol=1e-8)
        differences = (std_positions @ w.T + t)[:, None] - std_positions
        differences -= np.rint(differences)
        matches = (np.linalg.norm(differences @ std_basis.T, axis=-1) < 1e-8) & np.equal.outer(
            std_numbers, std_numbers
        )
        assert np.all(matches.sum(axis=1) == 1)

    for species in [1, 2]:
        # Equivalent origin choices can exchange a/b (and f/g in rutile).
        # Check each reported letter against its defining coordinates instead
        # of assuming that the algorithm chooses the input origin.
        letters = set(np.array(dataset.wyckoffs)[numbers == species])
        assert len(letters) == 1
        letter = letters.pop()
        assert set(np.array(dataset.site_symmetry_symbols)[numbers == species]) == {
            site_symmetries[species]
        }
        sites = std_positions[std_numbers == species]
        assert len(sites) == np.count_nonzero(numbers == species)
        if number == 225:
            assert letter in {"a", "b"}
            representatives = [np.full(3, 0 if letter == "a" else 0.5)]
        elif species == 1:
            assert letter in {"a", "b"}
            representatives = [np.array([0, 0, 0 if letter == "a" else 0.5])]
        else:
            assert letter in {"f", "g"}
            # 4f: (u,u,0); 4g: (u,-u,0), modulo lattice translations.
            constraints = np.column_stack(
                (sites[:, 2], sites[:, 0] - (1 if letter == "f" else -1) * sites[:, 1])
            )
            candidates = sites[np.all(np.abs(constraints - np.rint(constraints)) < 1e-8, axis=1)]
            assert len(candidates) > 0
            # Every site satisfying the defining coordinates must generate the orbit.
            representatives = candidates
        for representative in representatives:
            orbit = rotations @ representative + translations
            differences = orbit[:, None] - sites
            differences -= np.rint(differences)
            matches = np.linalg.norm(differences, axis=-1) < 1e-8
            assert np.all(matches.sum(axis=1) == 1)
            assert np.all(matches.sum(axis=0) == len(reference) // len(sites))
