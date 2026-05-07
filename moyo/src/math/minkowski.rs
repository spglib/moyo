use itertools::Itertools;
use nalgebra::allocator::Allocator;
use nalgebra::{
    DMatrix, DVector, DefaultAllocator, Dim, DimName, Matrix2, Matrix3, OMatrix, OVector, U3,
    Vector3,
};

use super::cycle_checker::CycleChecker;
use super::elementary::swapping_column_matrix;

const EPS: f64 = 1e-8;

/// basis is column-wise
/// If the basis is already Minkowski reduced, return the input basis and identity matrix
pub fn minkowski_reduce(basis: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<i32>) {
    if is_minkowski_reduced(basis) {
        return (*basis, Matrix3::<i32>::identity());
    }

    let mut reduced_basis = *basis;
    let mut trans_mat = Matrix3::<i32>::identity();
    minkowski_reduce_greedy(U3, &mut reduced_basis, &mut trans_mat, 3);

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
fn minkowski_reduce_greedy<N: Dim + DimName>(
    dim: N,
    basis: &mut OMatrix<f64, N, N>,
    trans_mat: &mut OMatrix<i32, N, N>,
    rank: usize,
) where
    DefaultAllocator: Allocator<N, N> + Allocator<N, N> + Allocator<N>,
{
    // Line 1: exit condition
    if rank == 1 {
        return;
    }

    let mut cc = CycleChecker::new();
    loop {
        // Line 3: sort basis vectors by their lengths
        let lengths: Vec<f64> = basis.column_iter().map(|column| column.norm()).collect();
        for i in 0..rank {
            for j in 0..(rank - 1 - i) {
                if lengths[j] > lengths[j + 1] + EPS {
                    basis.swap_columns(j, j + 1);
                    *trans_mat *= swapping_column_matrix(dim, j, j + 1);
                }
            }
        }

        // Line 4: Recursive call
        minkowski_reduce_greedy(dim, basis, trans_mat, rank - 1);

        // Line 5: Solve the closest vector problem (CVP) for basis.column(d - 1)
        // linear projection (Gram-Schmidt)
        let h = DMatrix::from_fn(rank - 1, rank - 1, |i, j| {
            basis.column(i).dot(&basis.column(j)) / basis.column(i).norm_squared()
        });
        let u = DVector::from_fn(rank - 1, |i, _| {
            basis.column(i).dot(&basis.column(rank - 1)) / basis.column(i).norm_squared()
        });
        let gs_coeffs = h.try_inverse().unwrap() * u;
        let gs_coeffs_rint = DVector::from_fn(rank - 1, |i, _| gs_coeffs[(i, 0)].round() as i32);

        // Since basis[0..(rank-1)] are already Minkowski reduced, we only need to check Voronoi relevant vectors from y_int
        let mut cvp_min = f64::INFINITY;
        let mut coeffs_argmin = DVector::<i32>::zeros(rank - 1);
        let mut c_argmin = OVector::<f64, N>::zeros();

        // For up to three dimensions, `[-1, 0, 1]` is sufficient once
        // `basis[0..(rank-1)]` is already Minkowski reduced, so this bounded
        // integer search is exact under that precondition.
        for voronoi_vector in (0..rank - 1).map(|_| -1..=1).multi_cartesian_product() {
            let coeffs = gs_coeffs_rint.clone() + DVector::from_iterator(rank - 1, voronoi_vector);
            let c = basis.columns(0, rank - 1) * coeffs.map(|e| e as f64);
            let cvp = (c.clone() - basis.column(rank - 1)).norm();
            if cvp < cvp_min {
                cvp_min = cvp;
                coeffs_argmin = coeffs.clone();
                c_argmin = c;
            }
        }

        // Line 6: update basis.column(rank - 1)
        for j in 0..3 {
            basis[(j, rank - 1)] -= c_argmin[j];
        }
        let mut add_mat = OMatrix::<i32, N, N>::identity();
        for i in 0..(rank - 1) {
            add_mat[(i, rank - 1)] = -coeffs_argmin[i];
        }
        *trans_mat *= add_mat;

        // Line 7: loop until length ordering is changed
        if basis.column(rank - 1).norm() + EPS > basis.column(rank - 2).norm() {
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

// Lagrange-Gauss provably terminates in O(log(max norm)) steps; the bound is
// only a defensive guard against pathological f64 inputs.
const MINKOWSKI_2D_MAX_ITERATIONS: usize = 1024;

/// Lagrange-Gauss reduction of a 2D basis (column-vector convention).
/// Returns the reduced basis and the unimodular transformation `T` such that
/// `basis * T == reduced` (exact under integer linear combinations of input columns).
/// Parity is preserved (det(T) > 0). If the input is already reduced, returns
/// `(basis, identity)`.
#[allow(dead_code)]
pub fn minkowski_reduce_2d(basis: &Matrix2<f64>) -> (Matrix2<f64>, Matrix2<i32>) {
    if is_minkowski_reduced_2d(basis) {
        return (*basis, Matrix2::<i32>::identity());
    }

    let mut b = *basis;
    let mut t = Matrix2::<i32>::identity();

    for _ in 0..MINKOWSKI_2D_MAX_ITERATIONS {
        if b.column(0).norm() > b.column(1).norm() + EPS {
            b.swap_columns(0, 1);
            t = t * Matrix2::new(0, 1, 1, 0);
        }

        let denom = b.column(0).norm_squared();
        if denom < EPS {
            break;
        }
        let m = (b.column(0).dot(&b.column(1)) / denom).round() as i32;
        if m == 0 {
            break;
        }

        let new_b2 = b.column(1) - (m as f64) * b.column(0);
        b.set_column(1, &new_b2);
        t = t * Matrix2::new(1, -m, 0, 1);
    }

    // Preserve parity: det(T) > 0. In 2D, negating both columns leaves det unchanged,
    // so flip a single column instead (negate b2 / second column of T).
    if t.map(|e| e as f64).determinant() < 0. {
        let new_b2 = -b.column(1);
        b.set_column(1, &new_b2);
        t = t * Matrix2::new(1, 0, 0, -1);
    }

    (b, t)
}

/// Lift a 2x2 unimodular transformation (acting on the in-plane axes) to a 3x3
/// transformation that leaves the third axis untouched. Used to compose the
/// in-plane reduction with the layer-group block form W_33 = +/- 1 (paper eq. 4).
#[allow(dead_code)]
pub fn lift_2d_to_3d(t: &Matrix2<i32>) -> Matrix3<i32> {
    Matrix3::new(t[(0, 0)], t[(0, 1)], 0, t[(1, 0)], t[(1, 1)], 0, 0, 0, 1)
}

/// Returns true iff the 2D basis is Minkowski reduced
/// (`|b1| <= |b2|` and `|2 b1 . b2| <= |b1|^2`).
#[allow(dead_code)]
pub fn is_minkowski_reduced_2d(basis: &Matrix2<f64>) -> bool {
    let n1 = basis.column(0).norm();
    let n2 = basis.column(1).norm();
    if n1 > n2 + EPS {
        return false;
    }
    let dot12 = basis.column(0).dot(&basis.column(1));
    if 2.0 * dot12.abs() > basis.column(0).norm_squared() + EPS {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix2, Matrix3, Vector2, Vector3};
    use rand::SeedableRng;
    use rand::prelude::*;
    use rand::rngs::StdRng;

    use super::{
        is_minkowski_reduced, is_minkowski_reduced_2d, lift_2d_to_3d, minkowski_reduce,
        minkowski_reduce_2d,
    };

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
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
            );
            let (reduced_basis, trans_mat) = minkowski_reduce(&basis);
            assert!(is_minkowski_reduced(&reduced_basis));
            assert_relative_eq!(basis * trans_mat.map(|e| e as f64), reduced_basis);
        }

        for _ in 0..256 {
            let basis = Matrix3::<f64>::new(
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
                rng.random(),
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

    #[test]
    fn test_minkowski_idempotent() {
        let basis = Matrix3::<f64>::from_columns(&[
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-0.5, 3.0_f64.sqrt() / 2.0, 0.0),
            Vector3::new(0.0, 0.0, 10.0),
        ]);
        assert!(is_minkowski_reduced(&basis));

        let (reduced_basis, _) = minkowski_reduce(&basis);
        assert_relative_eq!(reduced_basis, basis);
    }

    #[test]
    fn test_minkowski_2d_basic() {
        let basis = Matrix2::<f64>::from_columns(&[Vector2::new(1.0, 0.0), Vector2::new(0.0, 1.0)]);
        let (reduced, t) = minkowski_reduce_2d(&basis);
        assert!(is_minkowski_reduced_2d(&reduced));
        assert_relative_eq!(reduced, basis);
        assert_eq!(t, Matrix2::<i32>::identity());

        let basis = Matrix2::<f64>::from_columns(&[Vector2::new(1.0, 0.0), Vector2::new(5.0, 1.0)]);
        let (reduced, t) = minkowski_reduce_2d(&basis);
        assert!(is_minkowski_reduced_2d(&reduced));
        assert_relative_eq!(reduced, basis * t.map(|e| e as f64));
    }

    #[test]
    fn test_minkowski_2d_random() {
        let mut rng: StdRng = SeedableRng::from_seed([1; 32]);

        for _ in 0..512 {
            let basis = Matrix2::<f64>::new(
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
                rng.random::<i8>() as f64,
            );
            // Skip degenerate (rank-deficient) cases where reduction is undefined.
            if basis.determinant().abs() < 1e-6 {
                continue;
            }
            let (reduced, t) = minkowski_reduce_2d(&basis);
            assert!(is_minkowski_reduced_2d(&reduced));
            assert_relative_eq!(basis * t.map(|e| e as f64), reduced);
            // Idempotence.
            let (reduced2, t2) = minkowski_reduce_2d(&reduced);
            assert_relative_eq!(reduced2, reduced);
            assert_eq!(t2, Matrix2::<i32>::identity());
            // Parity.
            assert!(t.map(|e| e as f64).determinant() > 0.0);
        }
    }

    #[test]
    fn test_minkowski_2d_swap_invariance() {
        let basis = Matrix2::<f64>::from_columns(&[Vector2::new(3.0, 0.0), Vector2::new(1.0, 2.0)]);
        let (r1, _) = minkowski_reduce_2d(&basis);

        let mut swapped = basis;
        swapped.swap_columns(0, 1);
        let (r2, _) = minkowski_reduce_2d(&swapped);
        // Reduced bases agree up to column swap and column sign flip.
        let mut equal_pair = false;
        for s0 in [-1.0_f64, 1.0] {
            for s1 in [-1.0_f64, 1.0] {
                let cand = Matrix2::from_columns(&[s0 * r2.column(0), s1 * r2.column(1)]);
                if (cand - r1).norm() < 1e-6 {
                    equal_pair = true;
                }
                let cand_swap = Matrix2::from_columns(&[s0 * r2.column(1), s1 * r2.column(0)]);
                if (cand_swap - r1).norm() < 1e-6 {
                    equal_pair = true;
                }
            }
        }
        assert!(equal_pair);
    }

    #[test]
    fn test_lift_2d_to_3d() {
        let t = Matrix2::<i32>::new(2, -1, 1, 0);
        let lifted = lift_2d_to_3d(&t);
        assert_eq!(lifted, Matrix3::<i32>::new(2, -1, 0, 1, 0, 0, 0, 0, 1));
    }
}
