use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, OMatrix};

use super::elementary::{
    adding_column_matrix, adding_row_matrix, changing_column_sign_matrix, swapping_column_matrix,
    swapping_row_matrix,
};

/// Hermite normal form of MxN matrix such that h = l * basis * r
#[derive(Debug)]
pub struct SNF<M: Dim, N: Dim>
where
    DefaultAllocator: Allocator<i32, M, N> + Allocator<i32, M, M> + Allocator<i32, N, N>,
{
    pub h: OMatrix<i32, M, N>,
    pub l: OMatrix<i32, M, M>,
    pub r: OMatrix<i32, N, N>,
}

pub fn snf<M: Dim, N: Dim>(basis: &OMatrix<i32, M, N>) -> SNF<M, N>
where
    DefaultAllocator: Allocator<i32, M, N> + Allocator<i32, M, M> + Allocator<i32, N, N>,
{
    let (m, n) = basis.shape_generic();
    let mut h = basis.clone();
    let mut l = OMatrix::identity_generic(m, m);
    let mut r = OMatrix::identity_generic(n, n);

    // Process the `s`th column and row
    for s in 0..usize::min(m.value(), n.value()) {
        if (s..m.value()).all(|i| h[(i, s)] == 0) && (s..n.value()).all(|j| h[(s, j)] == 0) {
            continue;
        }

        loop {
            // Choose pivot element with the smallest absolute value
            let pivot = (s..m.value())
                .flat_map(|i| (s..n.value()).map(move |j| (i, j)))
                .filter(|&(i, j)| h[(i, j)] != 0)
                .min_by_key(|&(i, j)| h[(i, j)].abs())
                .unwrap();

            // Move pivot element to (s, s)
            h.swap_rows(s, pivot.0);
            l = swapping_row_matrix(m, s, pivot.0) * l;
            h.swap_columns(s, pivot.1);
            r = r * swapping_column_matrix(n, s, pivot.1);

            // Guarantee that h[(s, s)] is positive
            if h[(s, s)] < 0 {
                for i in 0..m.value() {
                    h[(i, s)] *= -1;
                }
                r = r * changing_column_sign_matrix(n, s);
            }
            assert_ne!(h[(s, s)], 0);

            // Eliminate (*, s) entries
            let mut update = false;
            for i in 0..m.value() {
                if i == s {
                    continue;
                }
                let k = h[(i, s)] / h[(s, s)];

                if k != 0 {
                    update = true;
                    // h[(i, :)] -= k * h[(s, :)]
                    for j in 0..n.value() {
                        h[(i, j)] -= k * h[(s, j)];
                    }
                    l = adding_row_matrix(m, s, i, -k) * l;
                }
            }

            // Eliminate (s, *) entries
            for j in 0..n.value() {
                if j == s {
                    continue;
                }
                let k = h[(s, j)] / h[(s, s)];

                if k != 0 {
                    update = true;
                    // h[(:, j)] -= k * h[(:, s)]
                    for i in 0..m.value() {
                        h[(i, j)] -= k * h[(i, s)];
                    }
                    r = r * adding_column_matrix(n, s, j, -k);
                }
            }

            // Continue until updating
            if !update {
                break;
            }
        }
    }

    assert_eq!(h, l.clone() * basis * r.clone());
    SNF::<M, N> { h, l, r }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, SMatrix};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::snf;

    #[test]
    fn test_snf_small() {
        {
            let m = matrix![
                2, 4, 4;
                -6, 6, 12;
                10, -4, -16;
            ];
            let snf = snf(&m);
            assert_eq!(snf.h, matrix![2, 0, 0; 0, 6, 0; 0, 0, 12]);
        }
    }

    #[test]
    fn test_snf_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let m = SMatrix::<i32, 3, 3>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = snf(&m);

            let m = SMatrix::<i32, 5, 7>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = snf(&m);

            let m = SMatrix::<i32, 7, 5>::from_fn(|_, _| rng.gen_range(-4..4));
            let _ = snf(&m);
        }
    }
}
