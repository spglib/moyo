use nalgebra::{Matrix2, Matrix3};
use serde::{Deserialize, Serialize};

use crate::math::{is_minkowski_reduced_2d, minkowski_reduce_2d};

use super::error::MoyoError;

/// Two-dimensional lattice (column-vector basis convention) used by the
/// layer-group pipeline for the in-plane block of a [`LayerLattice`].
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
    /// matrix.
    ///
    /// Caller must ensure the in-plane axes lie in the xy-plane (i.e. the
    /// third row of `basis_3d` is zero). [`crate::base::LayerLattice::new`]
    /// validates this; a `debug_assert!` here catches direct callers that
    /// forget. A non-zero third row would silently drop information --
    /// see PR #304 review comment 3198817773.
    pub fn from_inplane_of(basis_3d: &Matrix3<f64>) -> Self {
        debug_assert!(
            basis_3d[(2, 0)].abs() < 1e-6 && basis_3d[(2, 1)].abs() < 1e-6,
            "in-plane axes must be in the xy-plane (a_z = b_z = 0); got a_z={}, b_z={}",
            basis_3d[(2, 0)],
            basis_3d[(2, 1)],
        );
        Self::from_basis(Matrix2::new(
            basis_3d[(0, 0)],
            basis_3d[(0, 1)],
            basis_3d[(1, 0)],
            basis_3d[(1, 1)],
        ))
    }

    /// Lagrange-Gauss-reduce this 2D lattice. Returns the reduced lattice and
    /// the unimodular transformation `T` such that `self.basis * T == reduced.basis`.
    /// Mirrors [`crate::base::Lattice::minkowski_reduce`]: validates the result
    /// and returns [`MoyoError::MinkowskiReductionError`] if reduction failed.
    pub fn minkowski_reduce(&self) -> Result<(Self, Matrix2<i32>), MoyoError> {
        let (reduced_basis, trans_mat) = minkowski_reduce_2d(&self.basis);
        let reduced = Self {
            basis: reduced_basis,
        };
        if !reduced.is_minkowski_reduced() {
            return Err(MoyoError::MinkowskiReductionError);
        }
        Ok((reduced, trans_mat))
    }

    /// True iff the basis is 2D-Minkowski-reduced.
    pub fn is_minkowski_reduced(&self) -> bool {
        is_minkowski_reduced_2d(&self.basis)
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix2, Matrix3, matrix};

    use super::Lattice2D;

    #[test]
    fn test_from_inplane_of_extracts_upper_left_block() {
        // Third column (aperiodic axis) can be arbitrary; third row must be
        // zero to satisfy the layer-pipeline contract enforced by `LayerLattice`.
        let basis_3d: Matrix3<f64> = matrix![
            1.0, 2.0, 99.0;
            3.0, 4.0, 99.0;
            0.0, 0.0, 99.0;
        ];
        let l = Lattice2D::from_inplane_of(&basis_3d);
        assert_eq!(l.basis, Matrix2::new(1.0, 2.0, 3.0, 4.0));
    }

    #[test]
    fn test_minkowski_reduce_already_reduced_is_identity() {
        let l = Lattice2D::from_basis(Matrix2::new(1.0, 0.0, 0.0, 1.0));
        let (reduced, t) = l.minkowski_reduce().unwrap();
        assert_eq!(reduced.basis, l.basis);
        assert_eq!(t, Matrix2::new(1, 0, 0, 1));
        assert!(reduced.is_minkowski_reduced());
    }

    #[test]
    fn test_minkowski_reduce_skewed_round_trip() {
        let l = Lattice2D::from_basis(Matrix2::new(1.0, 4.0, 0.0, 1.0));
        let (reduced, t) = l.minkowski_reduce().unwrap();
        assert!(reduced.is_minkowski_reduced());
        let t_f = t.map(|e| e as f64);
        let reconstructed = l.basis * t_f;
        assert!((reconstructed - reduced.basis).norm() < 1e-12);
    }
}
