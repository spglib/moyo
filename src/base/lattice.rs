use nalgebra::base::{Matrix3, Vector3};

use crate::math::delaunay::delaunay_reduce;
use crate::math::minkowski::{is_minkowski_reduced, minkowski_reduce};
use crate::math::niggli::{is_niggli_reduced, niggli_reduce};

use super::error::MoyoError;
use super::transformation::{Linear, UnimodularLinear};

pub type ColumnBasis = Matrix3<f64>;

#[derive(Debug, Clone)]
pub struct Lattice {
    /// basis.column(i) is the i-th basis vector
    pub basis: ColumnBasis,
}

impl Lattice {
    pub fn new(basis: ColumnBasis) -> Self {
        Self { basis }
    }

    pub fn transform(&self, trans_mat: &Linear) -> Self {
        Self {
            basis: self.basis * trans_mat,
        }
    }

    pub fn transform_unimodular(&self, trans_mat: &UnimodularLinear) -> Self {
        Self {
            basis: self.basis * trans_mat.map(|e| e as f64),
        }
    }

    pub fn minkowski_reduce(&self) -> Result<(Self, UnimodularLinear), MoyoError> {
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

    pub fn niggli_reduce(&self) -> Result<(Self, UnimodularLinear), MoyoError> {
        let (reduced_basis, trans_mat) = niggli_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };

        if !reduced_lattice.is_niggli_reduced() {
            return Err(MoyoError::NiggliReductionError);
        }

        Ok((reduced_lattice, trans_mat))
    }

    /// Return true if basis vectors are Niggli reduced
    pub fn is_niggli_reduced(&self) -> bool {
        is_niggli_reduced(&self.basis)
    }

    pub fn delaunay_reduce(&self) -> Result<(Self, UnimodularLinear), MoyoError> {
        let (reduced_basis, trans_mat) = delaunay_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };

        Ok((reduced_lattice, trans_mat))
    }

    pub fn metric_tensor(&self) -> Matrix3<f64> {
        self.basis.transpose() * self.basis
    }

    pub fn cartesian_coords(&self, fractional_coords: &Vector3<f64>) -> Vector3<f64> {
        self.basis * fractional_coords
    }

    pub fn volume(&self) -> f64 {
        self.basis.determinant().abs()
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
            1.0, 1.0, -1.0;
            1.0, 0.0, 0.0;
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
