use std::collections::BTreeMap;

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use union_find::{QuickFindUf, UnionByRank, UnionFind};

use super::lattice::Lattice;
use super::permutation::Permutation;

/// Fractional coordinates
pub type Position = Vector3<f64>;
/// Atomic number
pub type AtomicSpecie = i32;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Representing a crystal structure
pub struct Cell {
    /// Lattice of the cell.
    pub lattice: Lattice,
    /// `positions[i]` is a fractional coordinates of the i-th site.
    pub positions: Vec<Position>,
    /// `numbers[i]` is an atomic number of the i-th site.
    pub numbers: Vec<AtomicSpecie>,
}

impl Cell {
    pub fn new(lattice: Lattice, positions: Vec<Position>, numbers: Vec<AtomicSpecie>) -> Self {
        if positions.len() != numbers.len() {
            panic!("positions and numbers should be the same length");
        }
        Self {
            lattice,
            positions,
            numbers,
        }
    }

    /// Return the number of atoms in the cell.
    pub fn num_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Rotate the cell by the given rotation matrix.
    pub fn rotate(&self, rotation_matrix: &Matrix3<f64>) -> Self {
        Self::new(
            self.lattice.rotate(rotation_matrix),
            self.positions.clone(),
            self.numbers.clone(),
        )
    }
}

/// If and only if the `i`th and `j`th atoms are equivalent, `orbits[i] == orbits[j]`.
/// For each orbit, only one of them satisfies `orbits[i] == i`.
pub fn orbits_from_permutations(num_atoms: usize, permutations: &[Permutation]) -> Vec<usize> {
    let mut uf = QuickFindUf::<UnionByRank>::new(num_atoms);
    for permutation in permutations.iter() {
        for i in 0..num_atoms {
            uf.union(i, permutation.apply(i));
        }
    }
    let mut identifier_mapping = BTreeMap::new();
    for i in 0..num_atoms {
        identifier_mapping.entry(uf.find(i)).or_insert(i);
    }

    (0..num_atoms)
        .map(|i| *identifier_mapping.get(&uf.find(i)).unwrap())
        .collect()
}

#[cfg(test)]
mod tests {
    use std::panic;

    use nalgebra::{vector, Matrix3};

    use super::{orbits_from_permutations, Cell};
    use crate::base::lattice::Lattice;
    use crate::base::permutation::Permutation;

    #[test]
    fn test_orbits_from_permutations() {
        {
            let num_atoms = 3;
            let permutations = vec![Permutation::new(vec![2, 1, 0])];
            assert_eq!(
                orbits_from_permutations(num_atoms, &permutations),
                vec![0, 1, 0]
            );
        }
        {
            let num_atoms = 3;
            let permutations = vec![Permutation::new(vec![1, 0, 2])];
            assert_eq!(
                orbits_from_permutations(num_atoms, &permutations),
                vec![0, 0, 2]
            );
        }
    }

    #[test]
    fn test_mismatched_length() {
        let lattice = Lattice::new(Matrix3::<f64>::identity());
        let positions = vec![vector![0.0, 0.0, 0.0], vector![0.5, 0.5, 0.5]];
        let numbers = vec![1];

        let result = panic::catch_unwind(|| Cell::new(lattice, positions, numbers));
        assert!(result.is_err());
    }
}
