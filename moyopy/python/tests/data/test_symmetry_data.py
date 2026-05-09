from __future__ import annotations

import numpy as np

from moyopy import Operations, operations_from_layer_number, operations_from_number


def _unique_sites_in_cell(position, operations: Operations) -> np.ndarray:
    def _is_close(site1, site2):
        diff = site1 - site2
        diff -= np.round(diff)
        return np.allclose(diff, 0, atol=1e-4)

    sites = []
    for rotation, translation in zip(operations.rotations, operations.translations):
        new_position = np.array(rotation) @ np.array(position) + np.array(translation)
        if np.any([_is_close(site, new_position) for site in sites]):
            continue
        sites.append(new_position)
    return np.array(sites)


def test_operations_from_number():
    # Test with default primitive=False
    operations = operations_from_number(number=230)  # Ia-3d
    num_operations = 48 * 2
    assert operations.num_operations == num_operations
    assert len(operations) == num_operations

    assert len(_unique_sites_in_cell([1 / 8, 1 / 8, 1 / 8], operations)) == 16

    # Test with primitive=True
    prim_operations = operations_from_number(number=230, primitive=True)  # Ia-3d
    assert prim_operations.num_operations == 48
    assert len(prim_operations) == 48


def test_operations_from_layer_number_primitive():
    # LG 1 (p 1) -- single identity operation in the primitive cell.
    operations = operations_from_layer_number(number=1)
    assert operations.num_operations == 1

    # LG 80 (p 6/m m m) -- order 24, primitive == conventional (centering P).
    operations = operations_from_layer_number(number=80)
    assert operations.num_operations == 24
    prim_operations = operations_from_layer_number(number=80, primitive=True)
    assert prim_operations.num_operations == 24


def test_operations_from_layer_number_centered():
    # LG 10 (c 2 1 1) -- centered (oc) rectangular.
    # Conventional cell has 2 lattice points x order-2 point group = 4 ops.
    operations = operations_from_layer_number(number=10)
    assert operations.num_operations == 4
    # Primitive cell has the centering folded out -> 2 ops.
    prim_operations = operations_from_layer_number(number=10, primitive=True)
    assert prim_operations.num_operations == 2
