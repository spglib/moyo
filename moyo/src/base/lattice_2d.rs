use nalgebra::{Matrix2, Matrix3};
use serde::{Deserialize, Serialize};

use crate::math::{is_minkowski_reduced_2d, minkowski_reduce_2d};

/// Two-dimensional lattice (column-vector basis convention) used by the
/// layer-group pipeline for the in-plane block of a [`LayerLattice`].
///
/// Kept separate from the bulk [`Lattice`] type because the layer pipeline
/// operates on the in-plane (a, b) sublattice independently of the aperiodic
/// `c` axis -- mixing them would defeat the layer-vs-bulk type separation
/// established by the [`LayerLattice`] / [`LayerCell`] newtypes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lattice2D {
    /// `basis.column(i)` is the i-th in-plane basis vector.
    pub basis: Matrix2<f64>,
}

impl Lattice2D {
    /// Construct from a column-basis 2x2 matrix.
    pub fn from_basis(basis: Matrix2<f64>) -> Self {
        Self { basis }
    }

    /// Construct from the upper-left 2x2 in-plane block of a 3D column-basis
    /// matrix, dropping the third row and column. Layer cells satisfy
    /// `c . a = c . b = 0` (paper Fu et al. 2024 eq. 5), so the third row is
    /// `(0, 0, *)` and the third column is the aperiodic axis -- both
    /// irrelevant to in-plane lattice operations.
    pub fn from_inplane_of(basis_3d: &Matrix3<f64>) -> Self {
        Self::from_basis(Matrix2::new(
            basis_3d[(0, 0)],
            basis_3d[(0, 1)],
            basis_3d[(1, 0)],
            basis_3d[(1, 1)],
        ))
    }

    /// Lagrange-Gauss-reduce this 2D lattice. Returns the reduced lattice and
    /// the unimodular transformation `T` such that `self.basis * T == reduced.basis`.
    pub fn minkowski_reduce(&self) -> (Self, Matrix2<i32>) {
        let (reduced, trans_mat) = minkowski_reduce_2d(&self.basis);
        (Self { basis: reduced }, trans_mat)
    }

    /// True iff the basis is 2D-Minkowski-reduced.
    pub fn is_minkowski_reduced(&self) -> bool {
        is_minkowski_reduced_2d(&self.basis)
    }
}
