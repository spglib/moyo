use itertools::Itertools;
use nalgebra::base::{DMatrix, DVector, Matrix3, Vector3};

use crate::transformation::{OriginShift, Transformation, TransformationMatrix};

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

    pub fn minkowski_reduce(&self) -> (Self, Transformation) {
        let mut minkowski_basis = self.basis.clone();
        let mut trans_mat = TransformationMatrix::identity();
        minkowski_reduce_greedy(&mut minkowski_basis, &mut trans_mat, 3);
        let minkowski_lattice = Self {
            basis: minkowski_basis,
        };
        let trans = Transformation::new(trans_mat, OriginShift::zeros());

        assert_relative_eq!(
            self.transform(&trans).basis,
            minkowski_lattice.basis,
            epsilon = EPS
        );
        assert!(minkowski_lattice.is_minkowski_reduced());

        (minkowski_lattice, trans)
    }

    /// Return true if basis vectors are Minkowski reduced
    pub fn is_minkowski_reduced(&self) -> bool {
        let norms = self.basis.column_iter().map(|v| v.norm()).collect_vec();

        // Check ordering: norms[0] <= norms[1] <= norms[2]
        if norms[0] > norms[1] + EPS {
            return false;
        }
        if norms[1] > norms[2] + EPS {
            return false;
        }

        // Check for norms[1]
        for coeffs in [[1, -1, 0], [1, 1, 0]] {
            let v = self.basis * Vector3::new(coeffs[0] as f64, coeffs[1] as f64, coeffs[2] as f64);
            if v.norm() + EPS < norms[1] {
                return false;
            }
        }

        // Check for norms[2]
        for coeffs in [
            [1, 0, 1],
            [1, 0, -1],
            [0, 1, 1],
            [0, 1, -1],
            [1, -1, -1],
            [1, -1, 1],
            [1, 1, -1],
            [1, 1, 1],
        ] {
            let v = self.basis * Vector3::new(coeffs[0] as f64, coeffs[1] as f64, coeffs[2] as f64);
            if v.norm() + EPS < norms[2] {
                return false;
            }
        }

        true
    }

    pub fn metric_tensor(&self) -> Matrix3<f64> {
        self.basis.transpose() * self.basis
    }
}

/// Implement Fig.3 in [1]
/// basis * trans_mat is always preserved
///
/// [1] https://royalsocietypublishing.org/doi/abs/10.1098/rspa.1992.0004
fn minkowski_reduce_greedy(
    basis: &mut ColumnBasis,
    trans_mat: &mut TransformationMatrix,
    dim: usize,
) {
    // Line 1: exit condition
    if dim == 1 {
        return;
    }

    loop {
        // Line 3: sort basis vectors by their lengths
        let lengths: Vec<f64> = basis.column_iter().map(|column| column.norm()).collect();
        for i in 0..dim {
            for j in 0..(dim - 1 - i) {
                if lengths[j] > lengths[j + 1] + EPS {
                    basis.swap_columns(j, j + 1);
                    *trans_mat = *trans_mat * swapping_column_matrix(j, j + 1);
                }
            }
        }

        // Line 4: Recursive call
        minkowski_reduce_greedy(basis, trans_mat, dim - 1);

        // Line 5: Solve the closest vector problem (CVP) for basis.column(d - 1)
        // linear projection (Gram-Schmidt)
        let h = DMatrix::from_fn(dim - 1, dim - 1, |i, j| {
            basis.column(i).dot(&basis.column(j)) / basis.column(i).norm_squared()
        });
        let u = DVector::from_fn(dim - 1, |i, _| {
            basis.column(i).dot(&basis.column(dim - 1)) / basis.column(i).norm_squared()
        });
        let gs_coeffs = h.try_inverse().unwrap() * u;
        let gs_coeffs_rint = DVector::from_fn(dim - 1, |i, _| gs_coeffs[(i, 0)].round() as i32);

        // Since basis[0..(dim-1)] are already Minkowski reduced, we only need to check Voronoi relevant vectors from y_int
        let mut cvp_min = f64::INFINITY;
        let mut coeffs_argmin = DVector::<i32>::zeros(dim - 1);
        let mut c_argmin = Vector3::<f64>::zeros();

        // [-1, 0, 1] is sufficient for up to three dimensional Minkowski-reduced basis
        for voronoi_vector in (0..dim - 1).map(|_| -1..=1).multi_cartesian_product() {
            let coeffs = gs_coeffs_rint.clone() + DVector::from_iterator(dim - 1, voronoi_vector);
            let c = basis.columns(0, dim - 1) * coeffs.map(|e| e as f64);
            let cvp = (c - basis.column(dim - 1)).norm();
            if cvp < cvp_min {
                cvp_min = cvp;
                coeffs_argmin = coeffs.clone();
                c_argmin = c.clone();
            }
        }

        // Line 6: update basis.column(dim - 1)
        for j in 0..3 {
            basis[(j, dim - 1)] -= c_argmin[j];
        }
        let mut add_mat = TransformationMatrix::identity();
        for i in 0..(dim - 1) {
            add_mat[(i, dim - 1)] = -coeffs_argmin[i];
        }
        *trans_mat = *trans_mat * add_mat;

        // Line 7: loop until length ordering is changed
        if basis.column(dim - 1).norm() + EPS > basis.column(dim - 2).norm() {
            break;
        }
    }
}

/// Return elementary matrix swapping the `col1`th and `col2`th columns
fn swapping_column_matrix(col1: usize, col2: usize) -> TransformationMatrix {
    let mut trans_mat = TransformationMatrix::zeros();
    for i in 0..3 {
        if i == col1 {
            trans_mat[(col1, col2)] = 1;
        } else if i == col2 {
            trans_mat[(col2, col1)] = 1;
        } else {
            trans_mat[(i, i)] = 1
        }
    }
    trans_mat
}

/// Return elementary matrix adding the `k`-multiplied `col1`th column into the `col2`th column
#[allow(dead_code)]
fn adding_column_matrix(col1: usize, col2: usize, k: i32) -> TransformationMatrix {
    let mut trans_mat = TransformationMatrix::identity();
    for i in 0..3 {
        if i == col1 {
            trans_mat[(col1, col2)] = k;
        }
    }
    trans_mat
}

#[cfg(test)]
mod tests {
    use crate::lattice::{adding_column_matrix, TransformationMatrix};

    use super::{swapping_column_matrix, ColumnBasis, Lattice};
    use nalgebra::{Matrix3, Vector3};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_is_minkowski_reduced() {
        let lattice = Lattice::new(ColumnBasis::from_columns(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]));
        assert!(lattice.is_minkowski_reduced());

        let lattice = Lattice::new(ColumnBasis::from_columns(&[
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
        ]));
        assert!(!lattice.is_minkowski_reduced());
    }

    #[test]
    fn test_swapping_column_matrix() {
        let actual = swapping_column_matrix(0, 1);
        let expect = Matrix3::new(0, 1, 0, 1, 0, 0, 0, 0, 1);
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_adding_column_matrix() {
        let actual = adding_column_matrix(0, 2, -1);
        let expect = Matrix3::new(1, 0, -1, 0, 1, 0, 0, 0, 1);
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_minkowski_reduction_small() {
        let lattice = Lattice::new(ColumnBasis::from_columns(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]));
        let (reduced_lattice, trans) = lattice.minkowski_reduce();
        assert_relative_eq!(reduced_lattice.basis, lattice.basis);
        assert_eq!(trans.trans_mat, TransformationMatrix::identity());

        let lattice = Lattice::new(ColumnBasis::from_columns(&[
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
        ]));
        let (reduced_lattice, _) = lattice.minkowski_reduce();
        assert_relative_eq!(
            reduced_lattice.basis,
            ColumnBasis::from_columns(&[
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ])
        );

        let lattice = Lattice::new(ColumnBasis::from_columns(&[
            Vector3::new(-5.0, -10.0, 17.0),
            Vector3::new(17.0, 24.0, 12.0),
            Vector3::new(-127.0, 73.0, 5.0),
        ]));
        let _ = lattice.minkowski_reduce();
    }

    #[test]
    fn debug_test_minkowski_reduction() {}

    #[test]
    fn test_minkowski_reduction_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let lattice = Lattice::new(ColumnBasis::new(
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
            ));
            let _ = lattice.minkowski_reduce();
        }

        for _ in 0..256 {
            let lattice = Lattice::new(ColumnBasis::new(
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
            ));
            let _ = lattice.minkowski_reduce();
        }
    }
}
