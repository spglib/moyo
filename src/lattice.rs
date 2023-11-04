use nalgebra::base::{Matrix3, Vector3};

use crate::error::MoyoError;
use crate::math::minkowski::{is_minkowski_reduced, minkowski_reduce};
use crate::transformation::{OriginShift, Transformation};

pub type ColumnBasis = Matrix3<f64>;

#[derive(Debug)]
pub struct Lattice {
    /// basis.column(i) is the i-th basis vector
    pub basis: ColumnBasis,
}

const EPS: f64 = 1e-8;

impl Lattice {
    pub fn new(basis: ColumnBasis) -> Self {
        Self { basis }
    }

    pub fn transform(&self, trans: &Transformation) -> Self {
        let trans_mat_as_f64 = trans.trans_mat_as_f64();
        Self {
            basis: self.basis * trans_mat_as_f64,
        }
    }

    pub fn minkowski_reduce(&self) -> Result<(Self, Transformation), MoyoError> {
        let (minkowski_basis, trans_mat) = minkowski_reduce(&self.basis);
        let minkowski_lattice = Self {
            basis: minkowski_basis,
        };
        let trans = Transformation::new(trans_mat, OriginShift::zeros());

        if relative_ne!(
            self.transform(&trans).basis,
            minkowski_lattice.basis,
            epsilon = EPS
        ) {
            return Err(MoyoError::MinkowskiReductionError);
        }
        if !minkowski_lattice.is_minkowski_reduced() {
            return Err(MoyoError::MinkowskiReductionError);
        }

        Ok((minkowski_lattice, trans))
    }

    /// Return true if basis vectors are Minkowski reduced
    pub fn is_minkowski_reduced(&self) -> bool {
        is_minkowski_reduced(&self.basis)
    }

    pub fn metric_tensor(&self) -> Matrix3<f64> {
        self.basis.transpose() * self.basis
    }

    pub fn cartesian_coords(&self, fractional_coords: &Vector3<f64>) -> Vector3<f64> {
        self.basis * fractional_coords
    }
}
