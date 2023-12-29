use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, DimMin, Dyn, OMatrix, OVector, U1};

use super::snf::SNF;

#[derive(Debug)]
pub struct IntegerLinearSystem<N: Dim>
where
    DefaultAllocator: Allocator<i32, N>,
{
    /// rank of the integer linear system
    pub rank: usize,
    /// Special solution for a * x = b
    pub x: OVector<i32, N>,
    /// Nullspace of a * x = 0
    pub nullspace: OMatrix<i32, Dyn, N>,
}

impl<N: Dim> IntegerLinearSystem<N>
where
    DefaultAllocator: Allocator<i32, N>,
{
    /// Solve a * x = b
    /// If no solution, return None
    pub fn new<M: DimMin<N>>(a: &OMatrix<i32, M, N>, b: &OVector<i32, M>) -> Option<Self>
    where
        DefaultAllocator:
            Allocator<i32, M, N> + Allocator<i32, M, M> + Allocator<i32, N, N> + Allocator<i32, M>,
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
