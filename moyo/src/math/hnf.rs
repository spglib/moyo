use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, OMatrix};

/// Hermite normal form of (M, N) matrix such that h = basis * r
#[derive(Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct HNF<M: Dim, N: Dim>
where
    DefaultAllocator: Allocator<M, N> + Allocator<N, N>,
{
    pub h: OMatrix<i32, M, N>,
    pub r: OMatrix<i32, N, N>,
}

impl<M: Dim, N: Dim> HNF<M, N>
where
    DefaultAllocator: Allocator<M, N> + Allocator<N, N>,
{
    /// Return column-wise Hermite norm form
    pub fn new(basis: &OMatrix<i32, M, N>) -> Self {
        let (m, n) = basis.shape_generic();
        let mut h = basis.clone();
        let mut r = OMatrix::identity_generic(n, n);

        // Process the `s`th row
        for s in 0..m.value() {
            loop {
                if (s..n.value()).all(|j| h[(s, j)] == 0) {
                    break;
                }

                // Choose pivot column with the smallest absolute value
                let pivot = (s..n.value())
                    .filter(|&j| h[(s, j)] != 0)
                    .min_by_key(|&j| h[(s, j)].abs())
                    .unwrap();
                h.swap_columns(s, pivot);
                r.swap_columns(s, pivot);

                // Guarantee that h[(s, s)] is positive
                if h[(s, s)] < 0 {
                    for i in 0..m.value() {
                        h[(i, s)] *= -1;
                    }
                    for i in 0..n.value() {
                        r[(i, s)] *= -1;
                    }
                }
                assert_ne!(h[(s, s)], 0);

                // Add the `s`th column to the other columns
                let mut update = false;
                for j in 0..n.value() {
                    if j == s {
                        continue;
                    }
                    let k = h[(s, j)].div_euclid(h[(s, s)]);

                    if k != 0 {
                        update = true;
                        // h[(:, j)] -= k * h[(:, s)]
                        for i in 0..m.value() {
                            h[(i, j)] -= k * h[(i, s)];
                        }
                        // r[(:, j)] -= k * r[(:, s)]
                        for i in 0..n.value() {
                            r[(i, j)] -= k * r[(i, s)];
                        }
                    }
                }

                // Continue until updating
                if !update {
                    break;
                }
            }
        }
        assert_eq!(h, basis * r.clone());
        Self { h, r }
    }
}

#[cfg(test)]
mod tests {
    use itertools::iproduct;
    use nalgebra::{Dyn, Matrix3, OMatrix, SMatrix, U3, Vector3, matrix};
    use rand::SeedableRng;
    use rand::prelude::*;
    use rand::rngs::StdRng;

    use super::HNF;

    #[test]
    fn test_hnf_small() {
        {
            let m = matrix![
                -1, 0, 0;
                1, 2, 2;
                0, -1, -2;
            ];
            let hnf = HNF::new(&m);
            let expect = matrix![
                1, 0, 0;
                1, 2, 0;
                0, 0, 1;
            ];
            assert_eq!(hnf.h, expect);
        }
        {
            let m = matrix![
                20, -6;
                -2, 1;
            ];
            let hnf = HNF::new(&m);
            assert_eq!(hnf.h, matrix![2, 0; 1, 4]);
        }
        {
            let m = matrix![
                2, 3, 6, 2;
                5, 6, 1, 6;
                8, 3, 1, 1;
            ];
            let hnf = HNF::new(&m);
            let expect = matrix![
                1, 0, 0, 0;
                0, 1, 0, 0;
                0, 0, 1, 0;
            ];
            assert_eq!(hnf.h, expect);
        }
    }

    #[test]
    fn test_hnf_wide() {
        // The 3 x (3 + n) shape `transformation_matrix_from_translations` builds, one
        // column per pure translation -- here the 256 of an fcc 4x4x4 supercell, scaled
        // by 256. Accumulating `r` by matrix product made this shape cost O(n^4).
        let n = 256;
        let mut columns = vec![
            Vector3::new(n, 0, 0),
            Vector3::new(0, n, 0),
            Vector3::new(0, 0, n),
        ];
        for (i, j, k) in iproduct!(0..4, 0..4, 0..4) {
            for (di, dj, dk) in [(0, 0, 0), (0, 32, 32), (32, 0, 32), (32, 32, 0)] {
                columns.push(Vector3::new(64 * i + di, 64 * j + dj, 64 * k + dk));
            }
        }
        let basis = OMatrix::<i32, U3, Dyn>::from_columns(&columns);
        assert_eq!(basis.ncols(), 3 + n as usize);

        let hnf = HNF::new(&basis);
        // Also asserted inside `HNF::new`; kept here in case that is ever relaxed.
        assert_eq!(hnf.h, basis * &hnf.r);

        // The leading block spans the translation lattice, of index 256.
        let leading = Matrix3::from_columns(&[hnf.h.column(0), hnf.h.column(1), hnf.h.column(2)]);
        assert_relative_eq!(
            leading.map(|e| e as f64).determinant(),
            (n as f64).powi(3) / n as f64
        );
        // The reduction clears everything past it.
        assert!((3..hnf.h.ncols()).all(|j| hnf.h.column(j).iter().all(|&e| e == 0)));
    }

    #[test]
    fn test_hnf_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let m = SMatrix::<i32, 3, 3>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = HNF::new(&m);

            let m = SMatrix::<i32, 5, 7>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = HNF::new(&m);

            let m = SMatrix::<i32, 7, 5>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = HNF::new(&m);
        }
    }
}
