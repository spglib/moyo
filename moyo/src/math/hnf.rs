use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, OMatrix};

use super::elementary::{
    adding_column_matrix, changing_column_sign_matrix, swapping_column_matrix,
};

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
                r *= swapping_column_matrix(n, s, pivot);

                // Guarantee that h[(s, s)] is positive
                if h[(s, s)] < 0 {
                    for i in 0..m.value() {
                        h[(i, s)] *= -1;
                    }
                    r *= changing_column_sign_matrix(n, s);
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
                        r *= adding_column_matrix(n, s, j, -k);
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
    use nalgebra::{matrix, SMatrix};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

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
