use std::collections::HashSet;

use nalgebra::base::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim, OMatrix};

/// Record transformation matrices during lattice reduction
#[derive(Debug)]
pub struct CycleChecker<N: Dim>
where
    DefaultAllocator: Allocator<N, N>,
{
    visited: HashSet<OMatrix<i32, N, N>>,
}

impl<N: Dim> CycleChecker<N>
where
    DefaultAllocator: Allocator<N, N>,
{
    pub fn new() -> Self {
        Self {
            visited: HashSet::new(),
        }
    }

    /// If `matrix` is not visited, insert it and return true.
    pub fn insert(&mut self, matrix: &OMatrix<i32, N, N>) -> bool {
        if self.visited.contains(matrix) {
            return false;
        }
        self.visited.insert(matrix.clone());
        true
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::CycleChecker;

    #[test]
    fn test_cycle_checker() {
        let mut cc = CycleChecker::new();
        assert_eq!(cc.insert(&matrix![1, 0, 0; 0, 1, 0; 0, 0, 1]), true);
        assert_eq!(cc.insert(&matrix![1, 0, 0; 0, 1, 0; 0, 0, 1]), false);
    }
}
