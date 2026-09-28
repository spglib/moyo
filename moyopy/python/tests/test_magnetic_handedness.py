from __future__ import annotations

import numpy as np
import pytest
from _helpers import REBASES, assert_proper_rotation

from moyopy import (
    CollinearMagneticCell,
    MagneticSpaceGroupType,
    MoyoCollinearMagneticDataset,
    MoyoNonCollinearMagneticDataset,
    NonCollinearMagneticCell,
    magnetic_operations_from_uni_number,
)


def _operation_keys(operations, linear=None, origin=None, index=1):
    rotations = np.tile(operations.rotations, (index, 1, 1))
    translations = np.concatenate(
        [np.array(operations.translations) + [0, 0, k] for k in range(index)]
    )
    if linear is not None:
        inverse = np.linalg.inv(linear)
        translations = (rotations @ origin + translations - origin) @ inverse.T
        rotations = inverse @ rotations @ linear
        np.testing.assert_allclose(rotations, np.rint(rotations), atol=1e-8)
    return {
        (
            tuple(np.rint(rotation).astype(int).flatten()),
            tuple(np.rint(translation * 1e8).astype(int) % 10**8),
            time_reversal,
        )
        for rotation, translation, time_reversal in zip(
            rotations, translations, np.tile(operations.time_reversals, index)
        )
    }


@pytest.mark.parametrize(
    "uni_number,bns_number,index",
    [
        (667, "76.7", 1),
        (668, "76.8", 1),
        (669, "76.9", 1),
        (670, "76.10", 1),
        (3, "1.3", 1),
        (23, "5.16", 1),
        (932, "113.272", 1),
        (670, "76.10", 2),
    ],
)
@pytest.mark.parametrize("noncollinear", [False, True])
@pytest.mark.parametrize("is_axial", [False, True])
@pytest.mark.parametrize("rotate_basis", [False, True])
@pytest.mark.parametrize("rebase", REBASES)
def test_magnetic_passive_rebasing(
    uni_number, bns_number, index, noncollinear, is_axial, rotate_basis, rebase
):
    if uni_number == 3:
        basis = np.array([[3.0, 0.2, 0.1], [0.0, 4.0, 0.3], [0.0, 0.0, 5.0]])
    elif uni_number == 23:
        basis = np.array([[3.0, 0.0, -1.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]])
    else:
        basis = np.diag([5.0, 5.0, 6.0])
    reference = magnetic_operations_from_uni_number(uni_number)
    positions, numbers, moments = [], [], []
    for species, seed, moment in [
        (1, [0.137, 0.239, 0.371], [0.2, 0.3, 1.0] if noncollinear else 1.0),
        (2, [0.291, 0.113, 0.083], [0.4, -0.2, 0.7] if noncollinear else 2.0),
    ]:
        # Type II has pure time reversal, so all moments must vanish.
        moment = np.asarray(moment) * (0 if uni_number == 668 else 1)
        for w, t, time_reversal in zip(
            reference.rotations, reference.translations, reference.time_reversals
        ):
            position = (np.array(w) @ seed + t) % 1
            if any(
                number == species
                and np.linalg.norm(position - other - np.rint(position - other)) < 1e-8
                for number, other in zip(numbers, positions)
            ):
                continue
            cartesian = basis @ w @ np.linalg.inv(basis)
            transformed_moment = cartesian @ moment if noncollinear else moment
            if is_axial:
                transformed_moment = transformed_moment * round(np.linalg.det(cartesian))
            if time_reversal:
                transformed_moment = -transformed_moment
            positions.append(position)
            numbers.append(species)
            moments.append(transformed_moment)
    positions, numbers, moments = np.array(positions), np.array(numbers), np.array(moments)
    # Build an explicit supercell using known translation cosets.
    supercell = np.diag([1, 1, index])
    positions = np.concatenate(
        [(positions + [0, 0, k]) @ np.linalg.inv(supercell).T for k in range(index)]
    )
    numbers = np.tile(numbers, index)
    moments = np.tile(moments, (index, 1)) if noncollinear else np.tile(moments, index)
    basis = basis @ supercell
    cartesian_rotation = np.array([[1.0, 0.0, 0.0], [0.0, 0.6, -0.8], [0.0, 0.8, 0.6]])
    basis = cartesian_rotation @ basis
    if noncollinear:
        moments = moments @ cartesian_rotation.T
    origin = np.array([0.13, 0.27, 0.19])
    input_basis = basis @ rebase
    input_positions = (positions - origin) @ np.linalg.inv(rebase).T
    np.testing.assert_allclose(
        input_positions @ input_basis.T + basis @ origin, positions @ basis.T, atol=1e-12
    )
    # Passive rebasing leaves the Cartesian moments untouched, even for axial moments.
    cell_type = NonCollinearMagneticCell if noncollinear else CollinearMagneticCell
    dataset_type = (
        MoyoNonCollinearMagneticDataset if noncollinear else MoyoCollinearMagneticDataset
    )
    cell = cell_type(
        input_basis.T.tolist(), input_positions.tolist(), numbers.tolist(), moments.tolist()
    )
    dataset = dataset_type(
        cell, symprec=1e-5, mag_symprec=1e-5, is_axial=is_axial, rotate_basis=rotate_basis
    )
    assert dataset.uni_number == uni_number
    assert MagneticSpaceGroupType(dataset.uni_number).bns_number == bns_number
    assert _operation_keys(dataset.magnetic_operations) == _operation_keys(
        reference, supercell @ rebase, supercell @ origin, index
    )
    np.testing.assert_array_equal(
        np.equal.outer(dataset.orbits, dataset.orbits), np.equal.outer(numbers, numbers)
    )
    rotation = np.array(dataset.std_rotation_matrix)
    assert_proper_rotation(rotation, rotate_basis)
    rotated_moments = moments @ rotation.T if noncollinear else moments
    for standardized, linear, shift, primitive in [
        (dataset.std_mag_cell, dataset.std_linear, dataset.std_origin_shift, False),
        (dataset.prim_std_mag_cell, dataset.prim_std_linear, dataset.prim_std_origin_shift, True),
    ]:
        assert np.linalg.det(standardized.basis) > 0
        assert np.sign(np.linalg.det(linear)) == np.sign(np.linalg.det(input_basis))
        np.testing.assert_allclose(
            np.array(standardized.basis).T, rotation @ input_basis @ linear, atol=1e-8
        )
        target = magnetic_operations_from_uni_number(uni_number, primitive=primitive)
        assert _operation_keys(
            dataset.magnetic_operations, linear, np.array(shift)
        ) == _operation_keys(target)
        expected_positions = (input_positions - shift) @ np.linalg.inv(linear).T
        diff = expected_positions[:, None] - np.array(standardized.positions)
        diff -= np.rint(diff)
        matches = (np.linalg.norm(diff, axis=-1) < 1e-7) & (
            numbers[:, None] == np.array(standardized.numbers)
        )
        assert np.all(matches.sum(axis=1) == 1)
        assert np.all(matches.sum(axis=0) == len(numbers) // len(standardized.numbers))
        for i, j in zip(*np.nonzero(matches)):
            np.testing.assert_allclose(
                standardized.magnetic_moments[j], rotated_moments[i], atol=1e-8
            )
        if primitive:
            np.testing.assert_array_equal(matches.argmax(axis=1), dataset.mapping_std_prim)
