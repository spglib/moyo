use itertools::Itertools;
use nalgebra::{DMatrix, DVector, Matrix3, Vector3, U3};

use super::cycle_checker::CycleChecker;
use super::elementary::swapping_column_matrix;

const EPS: f64 = 1e-8;

/// basis is column-wise
pub fn minkowski_reduce(basis: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<i32>) {
    let mut reduced_basis = *basis;
    let mut trans_mat = Matrix3::<i32>::identity();
    minkowski_reduce_greedy(&mut reduced_basis, &mut trans_mat, 3);

    // Preserve parity
    if trans_mat.map(|e| e as f64).determinant() < 0. {
        reduced_basis *= -1.;
        trans_mat *= -1;
    }

    (reduced_basis, trans_mat)
}

/// Implement Fig.3 in [1]
/// basis * trans_mat is always preserved
///
/// [1] https://royalsocietypublishing.org/doi/abs/10.1098/rspa.1992.0004
fn minkowski_reduce_greedy(basis: &mut Matrix3<f64>, trans_mat: &mut Matrix3<i32>, dim: usize) {
    // Line 1: exit condition
    if dim == 1 {
        return;
    }

    let mut cc = CycleChecker::new();
    loop {
        // Line 3: sort basis vectors by their lengths
        let lengths: Vec<f64> = basis.column_iter().map(|column| column.norm()).collect();
        for i in 0..dim {
            for j in 0..(dim - 1 - i) {
                if lengths[j] > lengths[j + 1] + EPS {
                    basis.swap_columns(j, j + 1);
                    *trans_mat *= swapping_column_matrix(U3, j, j + 1);
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
                c_argmin = c;
            }
        }

        // Line 6: update basis.column(dim - 1)
        for j in 0..3 {
            basis[(j, dim - 1)] -= c_argmin[j];
        }
        let mut add_mat = Matrix3::<i32>::identity();
        for i in 0..(dim - 1) {
            add_mat[(i, dim - 1)] = -coeffs_argmin[i];
        }
        *trans_mat *= add_mat;

        // Line 7: loop until length ordering is changed
        if basis.column(dim - 1).norm() + EPS > basis.column(dim - 2).norm() {
            break;
        }

        // If the new basis is already visited, stop the loop
        if !cc.insert(trans_mat) {
            break;
        }
    }
}

/// basis is column-wise
pub fn is_minkowski_reduced(basis: &Matrix3<f64>) -> bool {
    let norms = basis.column_iter().map(|v| v.norm()).collect_vec();

    // Check ordering: norms[0] <= norms[1] <= norms[2]
    if norms[0] > norms[1] + EPS {
        return false;
    }
    if norms[1] > norms[2] + EPS {
        return false;
    }

    // Check for norms[1]
    for coeffs in [[1, -1, 0], [1, 1, 0]] {
        let v = basis * Vector3::new(coeffs[0] as f64, coeffs[1] as f64, coeffs[2] as f64);
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
        let v = basis * Vector3::new(coeffs[0] as f64, coeffs[1] as f64, coeffs[2] as f64);
        if v.norm() + EPS < norms[2] {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Vector3};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::{is_minkowski_reduced, minkowski_reduce};

    #[test]
    fn test_is_minkowski_reduced() {
        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]);
        assert!(is_minkowski_reduced(&basis));

        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
        ]);
        assert!(!is_minkowski_reduced(&basis));
    }

    #[test]
    fn test_minkowski_reduction_small() {
        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ]);
        let (reduced_basis, trans_mat) = minkowski_reduce(&basis);
        assert_relative_eq!(reduced_basis, basis);
        assert_eq!(trans_mat, Matrix3::<i32>::identity());

        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 1.0),
        ]);
        let (reduced_basis, _) = minkowski_reduce(&basis);
        assert_relative_eq!(
            reduced_basis,
            Matrix3::<f64>::from_columns(&[
                Vector3::new(0.0, 1.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ])
        );

        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(-5.0, -10.0, 17.0),
            Vector3::new(17.0, 24.0, 12.0),
            Vector3::new(-127.0, 73.0, 5.0),
        ]);
        let _ = minkowski_reduce(&basis);
    }

    #[test]
    fn test_minkowski_reduction_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let basis = Matrix3::<f64>::new(
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
                rng.gen::<i8>() as f64,
            );
            let (reduced_basis, trans_mat) = minkowski_reduce(&basis);
            assert!(is_minkowski_reduced(&reduced_basis));
            assert_relative_eq!(basis * trans_mat.map(|e| e as f64), reduced_basis);
        }

        for _ in 0..256 {
            let basis = Matrix3::<f64>::new(
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
                rng.gen(),
            );
            let (reduced_basis, trans_mat) = minkowski_reduce(&basis);
            assert!(is_minkowski_reduced(&reduced_basis));
            assert_relative_eq!(
                basis * trans_mat.map(|e| e as f64),
                reduced_basis,
                epsilon = 1e-4
            );
        }
    }
}
