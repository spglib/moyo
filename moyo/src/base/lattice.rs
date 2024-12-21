use nalgebra::base::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::math::{
    delaunay_reduce, is_minkowski_reduced, is_niggli_reduced, minkowski_reduce, niggli_reduce,
};

use super::error::MoyoError;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Representing basis vectors of a lattice
pub struct Lattice {
    /// basis.column(i) is the i-th basis vector
    pub basis: Matrix3<f64>,
}

impl Lattice {
    /// Create a new lattice from row basis vectors
    pub fn new(row_basis: Matrix3<f64>) -> Self {
        Self {
            basis: row_basis.transpose(),
        }
    }

    /// Return Minkowski reduced lattice and transformation matrix to it
    pub fn minkowski_reduce(&self) -> Result<(Self, Matrix3<i32>), MoyoError> {
        let (reduced_basis, trans_mat) = minkowski_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };

        if !reduced_lattice.is_minkowski_reduced() {
            return Err(MoyoError::MinkowskiReductionError);
        }

        Ok((reduced_lattice, trans_mat))
    }

    /// Return true if basis vectors are Minkowski reduced
    pub fn is_minkowski_reduced(&self) -> bool {
        is_minkowski_reduced(&self.basis)
    }

    /// Return Niggli reduced lattice and transformation matrix to it
    pub fn niggli_reduce(&self) -> Result<(Self, Matrix3<i32>), MoyoError> {
        let (reduced_lattice, trans_mat) = self.unchecked_niggli_reduce();

        if !reduced_lattice.is_niggli_reduced() {
            return Err(MoyoError::NiggliReductionError);
        }

        Ok((reduced_lattice, trans_mat))
    }

    /// Return Niggli reduced lattice and transformation matrix to it without checking reduction condition
    pub fn unchecked_niggli_reduce(&self) -> (Self, Matrix3<i32>) {
        let (reduced_basis, trans_mat) = niggli_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };
        (reduced_lattice, trans_mat)
    }

    /// Return true if basis vectors are Niggli reduced
    pub fn is_niggli_reduced(&self) -> bool {
        is_niggli_reduced(&self.basis)
    }

    /// Return Delaunay reduced lattice and transformation matrix to it
    pub fn delaunay_reduce(&self) -> Result<(Self, Matrix3<i32>), MoyoError> {
        let (reduced_basis, trans_mat) = delaunay_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };

        Ok((reduced_lattice, trans_mat))
    }

    /// Return metric tensor of the basis vectors
    pub fn metric_tensor(&self) -> Matrix3<f64> {
        self.basis.transpose() * self.basis
    }

    /// Return cartesian coordinates from the given fractional coordinates
    pub fn cartesian_coords(&self, fractional_coords: &Vector3<f64>) -> Vector3<f64> {
        self.basis * fractional_coords
    }

    /// Return volume of the cell
    pub fn volume(&self) -> f64 {
        self.basis.determinant().abs()
    }

    /// Rotate the lattice by the given rotation matrix
    pub fn rotate(&self, rotation_matrix: &Matrix3<f64>) -> Self {
        Self {
            basis: rotation_matrix * self.basis,
        }
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::Lattice;

    #[test]
    fn test_metric_tensor() {
        let lattice = Lattice::new(matrix![
            1.0, 1.0, 1.0;
            1.0, 1.0, 0.0;
            1.0, -1.0, 0.0;
        ]);
        let metric_tensor = lattice.metric_tensor();
        assert_relative_eq!(
            metric_tensor,
            matrix![
                3.0, 2.0, 0.0;
                2.0, 2.0, 0.0;
                0.0, 0.0, 2.0;
            ]
        );
    }
}
