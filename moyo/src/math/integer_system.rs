use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, DimMin, Dyn, Matrix3, OMatrix, OVector, U1, U9};

use super::snf::SNF;

#[derive(Debug)]
pub struct IntegerLinearSystem<N: Dim>
where
    DefaultAllocator: Allocator<N>,
{
    /// rank of the integer linear system
    #[allow(dead_code)]
    pub rank: usize,
    /// Special solution for a * x = b
    #[allow(dead_code)]
    pub x: OVector<i32, N>,
    /// Nullspace of a * x = 0
    pub nullspace: OMatrix<i32, Dyn, N>,
}

impl<N: Dim> IntegerLinearSystem<N>
where
    DefaultAllocator: Allocator<N>,
{
    /// Solve a * x = b
    /// If no solution, return None
    pub fn new<M: DimMin<N>>(a: &OMatrix<i32, M, N>, b: &OVector<i32, M>) -> Option<Self>
    where
        DefaultAllocator: Allocator<M, N> + Allocator<M, M> + Allocator<N, N> + Allocator<M>,
    {
        let (_, n) = a.shape_generic();
        let snf = SNF::new(a);
        let rank = snf.rank();
        if rank == n.value() {
            return None;
        }

        // Solve snf.d * y = snf.l * b (x = snf.r * y)
        let mut y = OVector::zeros_generic(n, U1);
        let lb = snf.l * b;
        for i in 0..rank {
            if lb[(i, 0)] % snf.d[(i, i)] != 0 {
                return None;
            }
            y[i] = lb[(i, 0)] / snf.d[(i, i)];
        }
        let x = snf.r.clone() * y;

        let nullspace = snf.r.columns(rank, n.value() - rank).clone().transpose();

        Some(Self { rank, x, nullspace })
    }
}

/// Solve P^-1 * A[i] * P = B[i] (for all i)
/// vec(A * P - P * B) = (I_3 \otimes A - B^T \otimes I_3) * vec(P)
pub fn sylvester3(a: &[Matrix3<i32>], b: &[Matrix3<i32>]) -> Option<Vec<Matrix3<i32>>> {
    let size = a.len();
    assert_eq!(size, b.len());

    let mut coeffs = OMatrix::<i32, Dyn, U9>::zeros(9 * size);
    let identity = Matrix3::<i32>::identity();
    for k in 0..size {
        let adj = identity.kronecker(&a[k]) - b[k].transpose().kronecker(&identity);
        for i in 0..9 {
            for j in 0..9 {
                coeffs[(9 * k + i, j)] = adj[(i, j)];
            }
        }
    }
    let solution = IntegerLinearSystem::new(&coeffs, &OVector::<i32, Dyn>::zeros(coeffs.nrows()));

    if let Some(solution) = solution {
        let basis: Vec<_> = solution
            .nullspace
            .row_iter()
            .map(|e| {
                // Vectorization operator is column-major
                Matrix3::<i32>::new(
                    e[0], e[1], e[2], //
                    e[3], e[4], e[5], //
                    e[6], e[7], e[8], //
                )
                .transpose()
            })
            .collect();
        Some(basis)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};

    use super::IntegerLinearSystem;

    #[test]
    fn test_integer_linear_system() {
        {
            let a = matrix![
                6, 4, 10;
                -1, 1, -5;
            ];
            let b = vector![4, 11];
            let solution = IntegerLinearSystem::new(&a, &b).unwrap();
            assert_eq!(solution.rank, 2);
            assert_eq!(a * solution.x, b);
        }
        {
            let a = matrix![
                1, 1, 0;
            ];
            let b = vector![2];
            let solution = IntegerLinearSystem::new(&a, &b).unwrap();
            assert_eq!(solution.rank, 1);
            assert_eq!(a * solution.x, b);
        }
        {
            let a = matrix![
                2, 4, 0;
            ];
            let b = vector![1];
            let solution = IntegerLinearSystem::new(&a, &b);
            assert!(solution.is_none());
        }
    }
}
