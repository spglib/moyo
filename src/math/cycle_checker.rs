use std::collections::HashSet;

use nalgebra::Matrix3;

/// Record transformation matrices during lattice reduction
#[derive(Debug)]
pub struct CycleChecker {
    visited: HashSet<Matrix3<i32>>,
}

impl CycleChecker {
    pub fn new() -> Self {
        Self {
            visited: HashSet::new(),
        }
    }

    /// If `matrix` is not visited, insert it and return true.
    pub fn insert(&mut self, matrix: &Matrix3<i32>) -> bool {
        if self.visited.contains(matrix) {
            return false;
        }
        self.visited.insert(*matrix);
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
