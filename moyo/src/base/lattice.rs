use nalgebra::base::{Matrix3, OMatrix, RowVector3, Vector3};
use serde::{Deserialize, Serialize};

use crate::math::{
    delaunay_reduce, is_minkowski_reduced, is_niggli_reduced, minkowski_reduce, niggli_reduce,
};
use crate::utils::to_3x3_slice;

use super::error::MoyoError;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Representing basis vectors of a lattice
///
/// Three-dimensional reductions return a right-handed basis `B = self.basis * T`.
/// The integer transformation `T` is unimodular: its determinant is +1 for
/// right-handed input and -1 for left-handed input. Absolute volume is preserved.
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

    /// Create a new lattice from 3 row basis vectors
    pub fn from_basis(basis: [[f64; 3]; 3]) -> Self {
        Self::new(OMatrix::from_rows(&[
            RowVector3::from(basis[0]),
            RowVector3::from(basis[1]),
            RowVector3::from(basis[2]),
        ]))
    }

    /// Return a right-handed Minkowski reduced lattice and transformation matrix to it.
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

    /// Return true if basis vectors are Minkowski reduced, independently of handedness.
    pub fn is_minkowski_reduced(&self) -> bool {
        is_minkowski_reduced(&self.basis)
    }

    /// Return a right-handed Niggli reduced lattice and transformation matrix to it.
    pub fn niggli_reduce(&self) -> Result<(Self, Matrix3<i32>), MoyoError> {
        let (reduced_lattice, trans_mat) = self.unchecked_niggli_reduce();

        if !reduced_lattice.is_niggli_reduced() {
            return Err(MoyoError::NiggliReductionError);
        }

        Ok((reduced_lattice, trans_mat))
    }

    /// Return a right-handed Niggli reduced lattice and transformation matrix to it
    /// without checking the reduction condition.
    pub fn unchecked_niggli_reduce(&self) -> (Self, Matrix3<i32>) {
        let (reduced_basis, trans_mat) = niggli_reduce(&self.basis);
        let reduced_lattice = Self {
            basis: reduced_basis,
        };
        (reduced_lattice, trans_mat)
    }

    /// Return true if basis vectors are Niggli reduced, independently of handedness.
    pub fn is_niggli_reduced(&self) -> bool {
        is_niggli_reduced(&self.basis)
    }

    /// Return a right-handed Delaunay reduced lattice and transformation matrix to it.
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

    #[allow(dead_code)]
    pub(crate) fn lattice_constant(&self) -> [f64; 6] {
        let g = self.metric_tensor();
        let a = g[(0, 0)].sqrt();
        let b = g[(1, 1)].sqrt();
        let c = g[(2, 2)].sqrt();
        let alpha = (g[(1, 2)] / (b * c)).acos().to_degrees();
        let beta = (g[(0, 2)] / (a * c)).acos().to_degrees();
        let gamma = (g[(0, 1)] / (a * b)).acos().to_degrees();
        [a, b, c, alpha, beta, gamma]
    }

    /// Returns the basis vectors as a 3x3 array of row vectors.
    ///
    /// Each inner array `[x, y, z]` is a basis vector in Cartesian coordinates.
    /// `result[0]`, `result[1]`, `result[2]` correspond to the first, second,
    /// and third basis vectors respectively.
    ///
    /// This is a convenience method for users who do not depend on `nalgebra`.
    /// It is equivalent to transposing `self.basis` (which stores column vectors).
    pub fn basis_as_array(&self) -> [[f64; 3]; 3] {
        to_3x3_slice(&self.basis.transpose())
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
    use nalgebra::{Matrix3, matrix, vector};
    use rstest::rstest;

    use super::{Lattice, MoyoError};

    #[rstest]
    #[case::niggli(Lattice::niggli_reduce, true)]
    #[case::unchecked_niggli(|lattice: &Lattice| Ok(lattice.unchecked_niggli_reduce()), true)]
    #[case::delaunay(Lattice::delaunay_reduce, false)]
    #[case::minkowski(Lattice::minkowski_reduce, false)]
    fn test_reduction_handedness(
        #[case] reduce: fn(&Lattice) -> Result<(Lattice, Matrix3<i32>), MoyoError>,
        #[case] niggli: bool,
        #[values(1.0, -1.0)] handedness: f64,
    ) {
        for row_basis in [
            Matrix3::identity(),
            Matrix3::from_diagonal(&vector![1.0, 2.0, 10.0]),
            matrix![4.0, 0.0, 0.0; 3.0, 2.0, 0.0; 1.0, 1.0, 3.0],
        ] {
            let lattice =
                Lattice::new(Matrix3::from_diagonal(&vector![1.0, 1.0, handedness]) * row_basis);
            let (reduced, transformation) = reduce(&lattice).unwrap();
            assert_relative_eq!(
                reduced.basis,
                lattice.basis * transformation.map(f64::from),
                epsilon = 1e-12
            );
            let determinant = transformation
                .column(0)
                .cross(&transformation.column(1))
                .dot(&transformation.column(2));
            assert_eq!(determinant, handedness as i32);
            assert!(reduced.basis.determinant() > 0.0);
            assert_relative_eq!(reduced.volume(), lattice.volume(), epsilon = 1e-12);
            assert!(reduced.is_minkowski_reduced());
            if niggli {
                assert!(reduced.is_niggli_reduced());
            }
        }
    }

    #[test]
    fn test_left_handed_reduced_predicates_and_idempotence() {
        let lattice = Lattice::new(Matrix3::from_diagonal(&vector![1.0, 2.0, -3.0]));
        assert!(lattice.is_niggli_reduced());
        assert!(lattice.is_minkowski_reduced());
        for reduce in [Lattice::niggli_reduce, Lattice::minkowski_reduce] {
            let (reduced, transformation) = reduce(&lattice).unwrap();
            assert_relative_eq!(reduced.basis, -lattice.basis);
            assert_eq!(transformation, -Matrix3::<i32>::identity());
            let (twice_reduced, second_transformation) = reduce(&reduced).unwrap();
            assert_relative_eq!(twice_reduced.basis, reduced.basis);
            assert_eq!(second_transformation, Matrix3::<i32>::identity());
        }
    }

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
