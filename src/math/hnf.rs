use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, DimName, OMatrix};

use crate::math::elementary::{adding_column_matrix, changing_sign_matrix, swapping_column_matrix};

/// Hermite normal form of RxC matrix such that h = basis * r
#[derive(Debug)]
pub struct HNF<R: Dim, C: Dim>
where
    DefaultAllocator: Allocator<i32, R, C> + Allocator<i32, C, C>,
{
    pub h: OMatrix<i32, R, C>,
    pub r: OMatrix<i32, C, C>,
}

/// Return column-wise Hermite norm form
pub fn hnf<R: DimName, C: DimName>(basis: &OMatrix<i32, R, C>) -> HNF<R, C>
where
    DefaultAllocator: Allocator<i32, R, C> + Allocator<i32, C, C>,
{
    let mut h = basis.clone();
    let mut r = OMatrix::<i32, C, C>::identity();
    let m = R::dim();
    let n = C::dim();

    // Process the `s`th row
    for s in 0..m {
        if (s..n).all(|j| h[(s, j)] == 0) {
            continue;
        }

        loop {
            // Choose pivot column with the smallest absolute value
            let pivot = (s..n)
                .filter(|&j| h[(s, j)] != 0)
                .min_by_key(|&j| h[(s, j)].abs())
                .unwrap();
            h.swap_columns(s, pivot);
            r = r * swapping_column_matrix::<C>(s, pivot);

            // Guarantee that h[(s, s)] is positive
            if h[(s, s)] < 0 {
                for i in 0..m {
                    h[(i, s)] *= -1;
                }
                r = r * changing_sign_matrix(s);
            }
            assert_ne!(h[(s, s)], 0);

            // Add the `s`th column to the other columns
            let mut update = false;
            for j in 0..n {
                if j == s {
                    continue;
                }
                let k = h[(s, j)] / h[(s, s)];

                if k != 0 {
                    update = true;
                    // h[(:, j)] -= k * h[(:, s)]
                    for i in 0..m {
                        h[(i, j)] -= k * h[(i, s)];
                    }
                    r = r * adding_column_matrix(s, j, -k);
                }
            }

            // Continue until updating
            if !update {
                break;
            }
        }
    }
    assert_eq!(h, basis * r.clone());
    HNF::<R, C> { h, r }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, SMatrix};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::hnf;

    #[test]
    fn test_hnf_small() {
        {
            let m = matrix![
                20, -6;
                -2, 1;
            ];
            let hnf = hnf(&m);
            assert_eq!(hnf.h, matrix![2, 0; 1, 4]);
        }
        {
            let m = matrix![
                2, 3, 6, 2;
                5, 6, 1, 6;
                8, 3, 1, 1;
            ];
            let hnf = hnf(&m);
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
            let m = SMatrix::<i32, 3, 3>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = hnf(&m);

            let m = SMatrix::<i32, 5, 7>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = hnf(&m);

            let m = SMatrix::<i32, 7, 5>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = hnf(&m);
        }
    }
}
