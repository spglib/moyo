import numpy as np
import pytest

REBASES = [
    pytest.param(np.eye(3, dtype=int), id="identity"),
    pytest.param(np.diag([-1, 1, 1]), id="axis-flip"),
    pytest.param(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]), id="odd-permutation"),
    pytest.param(np.array([[0, -1, 2], [1, 0, 1], [0, 0, -1]]), id="skew-left"),
]


def assert_proper_rotation(rotation, rotate_basis):
    np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-8)
    assert np.linalg.det(rotation) == pytest.approx(1)
    if not rotate_basis:
        np.testing.assert_allclose(rotation, np.eye(3), atol=1e-8)
