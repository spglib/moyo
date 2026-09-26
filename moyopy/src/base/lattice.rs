use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::base::Lattice;
use moyo::utils::to_3x3_slice;

use super::PyMoyoError;

type LatticeReduction = ([[f64; 3]; 3], [[i32; 3]; 3]);

fn lattice_from_basis(basis: [[f64; 3]; 3]) -> PyResult<Lattice> {
    if !basis.iter().flatten().all(|value| value.is_finite()) {
        return Err(PyValueError::new_err(
            "basis must contain only finite values",
        ));
    }
    let lattice = Lattice::from_basis(basis);
    let determinant = lattice.basis.determinant();
    if !determinant.is_finite() || determinant == 0.0 {
        return Err(PyValueError::new_err(
            "basis must have a finite, nonzero determinant",
        ));
    }
    Ok(lattice)
}

/// Return a Niggli-reduced basis and its integer transformation matrix.
///
/// Parameters
/// ----------
/// basis : list\[list\[float\]\]
///     Three linearly independent, finite row-wise lattice vectors.
///
/// Returns
/// -------
/// reduced_basis : list\[list\[float\]\]
///     Row-wise basis vectors of the reduced lattice.
/// transformation : list\[list\[int\]\]
///     Unimodular matrix satisfying ``reduced_basis = transformation.T @ basis``
///     in NumPy notation.
///
/// Raises
/// ------
/// ValueError
///     If the basis is invalid or reduction fails.
#[pyfunction]
pub fn niggli_reduce(basis: [[f64; 3]; 3]) -> PyResult<LatticeReduction> {
    let (reduced, transformation) = lattice_from_basis(basis)?
        .niggli_reduce()
        .map_err(PyMoyoError::from)?;
    Ok((reduced.basis_as_array(), to_3x3_slice(&transformation)))
}

/// Return a Delaunay-reduced basis and its integer transformation matrix.
///
/// Parameters
/// ----------
/// basis : list\[list\[float\]\]
///     Three linearly independent, finite row-wise lattice vectors.
///
/// Returns
/// -------
/// reduced_basis : list\[list\[float\]\]
///     Row-wise basis vectors of the reduced lattice.
/// transformation : list\[list\[int\]\]
///     Unimodular matrix satisfying ``reduced_basis = transformation.T @ basis``
///     in NumPy notation.
///
/// Raises
/// ------
/// ValueError
///     If the basis is invalid or reduction fails.
#[pyfunction]
pub fn delaunay_reduce(basis: [[f64; 3]; 3]) -> PyResult<LatticeReduction> {
    let (reduced, transformation) = lattice_from_basis(basis)?
        .delaunay_reduce()
        .map_err(PyMoyoError::from)?;
    Ok((reduced.basis_as_array(), to_3x3_slice(&transformation)))
}

/// Return a Minkowski-reduced basis and its integer transformation matrix.
///
/// Parameters
/// ----------
/// basis : list\[list\[float\]\]
///     Three linearly independent, finite row-wise lattice vectors.
///
/// Returns
/// -------
/// reduced_basis : list\[list\[float\]\]
///     Row-wise basis vectors of the reduced lattice.
/// transformation : list\[list\[int\]\]
///     Unimodular matrix satisfying ``reduced_basis = transformation.T @ basis``
///     in NumPy notation.
///
/// Raises
/// ------
/// ValueError
///     If the basis is invalid or reduction fails.
#[pyfunction]
pub fn minkowski_reduce(basis: [[f64; 3]; 3]) -> PyResult<LatticeReduction> {
    let (reduced, transformation) = lattice_from_basis(basis)?
        .minkowski_reduce()
        .map_err(PyMoyoError::from)?;
    Ok((reduced.basis_as_array(), to_3x3_slice(&transformation)))
}

/// Return whether the row-wise basis vectors are Niggli reduced.
///
/// Raises ``ValueError`` if the basis is not finite and linearly independent.
#[pyfunction]
pub fn is_niggli_reduced(basis: [[f64; 3]; 3]) -> PyResult<bool> {
    Ok(lattice_from_basis(basis)?.is_niggli_reduced())
}

/// Return whether the row-wise basis vectors are Minkowski reduced.
///
/// Raises ``ValueError`` if the basis is not finite and linearly independent.
#[pyfunction]
pub fn is_minkowski_reduced(basis: [[f64; 3]; 3]) -> PyResult<bool> {
    Ok(lattice_from_basis(basis)?.is_minkowski_reduced())
}
