use nalgebra::base::{Matrix3, Vector3};

use crate::math::minkowski::{is_minkowski_reduced, minkowski_reduce};

use super::error::MoyoError;
use super::tolerance::EPS;
use super::transformation::TransformationMatrix;

pub type ColumnBasis = Matrix3<f64>;

#[derive(Debug, Clone)]
pub struct Lattice {
    /// basis.column(i) is the i-th basis vector
    pub basis: ColumnBasis,
}

#[allow(non_camel_case_types)]
enum BravaisTypeOfLattice {
    aP,
    mP,
    mC,
    oP,
    oS,
    oF,
    oI,
    tP,
    tI,
    hR,
    hP,
    cP,
    cF,
    cI,
}

impl Lattice {
    pub fn new(basis: ColumnBasis) -> Self {
        Self { basis }
    }

    pub fn transform(&self, trans_mat: &TransformationMatrix) -> Self {
        let trans_mat_as_f64 = trans_mat.map(|e| e as f64);
        Self {
            basis: self.basis * trans_mat_as_f64,
        }
    }

    pub fn minkowski_reduce(&self) -> Result<(Self, TransformationMatrix), MoyoError> {
        let (minkowski_basis, trans_mat) = minkowski_reduce(&self.basis);
        let minkowski_lattice = Self {
            basis: minkowski_basis,
        };

        if relative_ne!(
            self.transform(&trans_mat).basis,
            minkowski_lattice.basis,
            epsilon = EPS
        ) {
            return Err(MoyoError::MinkowskiReductionError);
        }
        if !minkowski_lattice.is_minkowski_reduced() {
            return Err(MoyoError::MinkowskiReductionError);
        }

        Ok((minkowski_lattice, trans_mat))
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
