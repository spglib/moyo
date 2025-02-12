use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, DimMin, OMatrix};

/// Hermite normal form of MxN matrix such that d = l * basis * r
#[derive(Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct SNF<M: DimMin<N>, N: Dim>
where
    DefaultAllocator: Allocator<M, N> + Allocator<M, M> + Allocator<N, N>,
{
    pub d: OMatrix<i32, M, N>,
    pub l: OMatrix<i32, M, M>,
    pub r: OMatrix<i32, N, N>,
}

impl<M: DimMin<N>, N: Dim> SNF<M, N>
where
    DefaultAllocator: Allocator<M, N> + Allocator<M, M> + Allocator<N, N>,
{
    pub fn new(basis: &OMatrix<i32, M, N>) -> SNF<M, N>
    where
        DefaultAllocator: Allocator<M, N> + Allocator<M, M> + Allocator<N, N>,
    {
        let (m, n) = basis.shape_generic();
        let mut d = basis.clone();
        let mut l = OMatrix::identity_generic(m, m);
        let mut r = OMatrix::identity_generic(n, n);

        // Process the `s`th column and row
        for s in 0..m.min(n).value() {
            // Choose pivot element with the nonzero smallest absolute value
            while let Some(pivot) = (s..m.value())
                .flat_map(|i| (s..n.value()).map(move |j| (i, j)))
                .filter(|&(i, j)| d[(i, j)] != 0)
                .min_by_key(|&(i, j)| d[(i, j)].abs())
            {
                // Move pivot element to (s, s)
                d.swap_rows(s, pivot.0);
                l.swap_rows(s, pivot.0);
                d.swap_columns(s, pivot.1);
                r.swap_columns(s, pivot.1);

                // Guarantee that h[(s, s)] is positive
                if d[(s, s)] < 0 {
                    for i in 0..m.value() {
                        d[(i, s)] *= -1;
                    }
                    // r *= changing_column_sign_matrix(n, s);
                    for i in 0..n.value() {
                        r[(i, s)] *= -1;
                    }
                }
                assert_ne!(d[(s, s)], 0);

                // Eliminate (*, s) entries
                let mut update = false;
                for i in s + 1..m.value() {
                    let k = d[(i, s)] / d[(s, s)];

                    if k != 0 {
                        update = true;
                        // d[(i, :)] -= k * d[(s, :)]
                        for j in 0..n.value() {
                            d[(i, j)] -= k * d[(s, j)];
                        }
                        // l[(i, :)] -= k * l[(s, :)]
                        for j in 0..m.value() {
                            l[(i, j)] -= k * l[(s, j)];
                        }
                    }
                }

                // Eliminate (s, *) entries
                for j in s + 1..n.value() {
                    let k = d[(s, j)] / d[(s, s)];

                    if k != 0 {
                        update = true;
                        // d[(:, j)] -= k * d[(:, s)]
                        for i in 0..m.value() {
                            d[(i, j)] -= k * d[(i, s)];
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

        assert_eq!(d, l.clone() * basis * r.clone());
        SNF::<M, N> { d, l, r }
    }

    pub fn rank(&self) -> usize {
        let (m, n) = self.d.shape_generic();
        (0..m.min(n).value())
            .filter(|&i| self.d[(i, i)] != 0)
            .count()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, SMatrix};
    use rand::prelude::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use super::SNF;

    #[test]
    fn test_snf_small() {
        {
            let m = matrix![
                2, 4, 4;
                -6, 6, 12;
                10, -4, -16;
            ];
            let snf = SNF::new(&m);
            assert_eq!(snf.d, matrix![2, 0, 0; 0, 6, 0; 0, 0, 12]);
        }
    }

    #[test]
    fn test_snf_random() {
        let mut rng: StdRng = SeedableRng::from_seed([0; 32]);

        for _ in 0..256 {
            let m = SMatrix::<i32, 3, 3>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = SNF::new(&m);

            let m = SMatrix::<i32, 5, 7>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = SNF::new(&m);

            let m = SMatrix::<i32, 7, 5>::from_fn(|_, _| rng.random_range(-4..4));
            let _ = SNF::new(&m);
        }
    }
}
