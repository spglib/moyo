from __future__ import annotations

import numpy as np

from moyopy import Operations, operations_from_number


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
    operations = operations_from_number(number=230)  # Ia-3d
    num_operations = 48 * 2
    assert operations.num_operations == num_operations
    assert len(operations) == num_operations

    assert len(_unique_sites_in_cell([1 / 8, 1 / 8, 1 / 8], operations)) == 16
