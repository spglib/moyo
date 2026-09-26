import numpy as np
import pytest

from moyopy import (
    delaunay_reduce,
    is_minkowski_reduced,
    is_niggli_reduced,
    minkowski_reduce,
    niggli_reduce,
)

REDUCTIONS = [niggli_reduce, delaunay_reduce, minkowski_reduce]
LATTICE_FUNCTIONS = [*REDUCTIONS, is_niggli_reduced, is_minkowski_reduced]


def test_niggli_triclinic_reference():
    # Reference reduced cell from spglib 2.7.0, away from reduction boundaries.
    basis = [[4.0, 0.0, 0.0], [3.0, 2.0, 0.0], [1.0, 1.0, 3.0]]
    reduced, _ = niggli_reduce(basis)
    np.testing.assert_allclose(reduced, [[-1.0, 2.0, 0.0], [1.0, 1.0, 3.0], [3.0, 2.0, 0.0]])


@pytest.mark.parametrize("reduce", REDUCTIONS)
@pytest.mark.parametrize(
    "basis",
    [
        pytest.param(np.eye(3), id="cubic"),
        pytest.param(np.diag([1.0, 2.0, 10.0]), id="elongated"),
        pytest.param([[4.0, 0.0, 0.0], [3.0, 2.0, 0.0], [1.0, 1.0, 3.0]], id="skew"),
        pytest.param([[4.0, 0.0, 0.0], [3.0, 2.0, 0.0], [1.0, 1.0, -3.0]], id="left-handed"),
    ],
)
def test_lattice_reduction(reduce, basis):
    basis = np.array(basis)
    input_basis = basis.tolist()
    reduced, transformation = reduce(basis=input_basis)

    assert isinstance(reduced, list)
    assert all(isinstance(value, float) for row in reduced for value in row)
    assert isinstance(transformation, list)
    assert all(isinstance(value, int) for row in transformation for value in row)
    reduced = np.array(reduced)
    transformation = np.array(transformation)
    np.testing.assert_allclose(reduced, transformation.T @ basis, atol=1e-12)
    assert round(np.linalg.det(transformation)) == 1
    np.testing.assert_allclose(np.linalg.det(reduced), np.linalg.det(basis))
    np.testing.assert_array_equal(input_basis, basis)

    # Every reduction produces a Minkowski-reduced basis in three dimensions.
    assert is_minkowski_reduced(reduced.tolist())
    if reduce is niggli_reduce:
        assert is_niggli_reduced(reduced.tolist())


@pytest.mark.parametrize(
    ("reduce", "is_reduced"),
    [(niggli_reduce, is_niggli_reduced), (minkowski_reduce, is_minkowski_reduced)],
)
def test_reduction_predicates_and_idempotence(reduce, is_reduced):
    basis = [[1.0, 0.0, 0.0], [5.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    assert not is_reduced(basis)
    reduced, _ = reduce(basis)
    assert is_reduced(reduced)
    twice_reduced, transformation = reduce(reduced)
    np.testing.assert_allclose(twice_reduced, reduced)
    np.testing.assert_array_equal(transformation, np.eye(3, dtype=int))


@pytest.mark.parametrize("function", LATTICE_FUNCTIONS)
@pytest.mark.parametrize(
    ("basis", "message"),
    [
        pytest.param(np.zeros((3, 3)), "nonzero determinant", id="zero"),
        pytest.param([[1, 0, 0], [2, 0, 0], [0, 0, 1]], "nonzero determinant", id="singular"),
        pytest.param(np.diag([1, 1, np.nan]), "finite values", id="nan"),
        pytest.param(np.diag([1, 1, np.inf]), "finite values", id="infinity"),
    ],
)
def test_invalid_lattice(function, basis, message):
    with pytest.raises(ValueError, match=message):
        function(np.asarray(basis).tolist())


@pytest.mark.parametrize("function", LATTICE_FUNCTIONS)
@pytest.mark.parametrize("basis", [[], [[1.0, 0.0], [0.0, 1.0]], [[1.0], [2.0], [3.0]]])
def test_invalid_basis_shape(function, basis):
    with pytest.raises(ValueError):
        function(basis)
