use nalgebra::Vector3;
use std::collections::BTreeMap;
use union_find::{QuickFindUf, UnionByRank, UnionFind};

use super::lattice::Lattice;
use super::operation::Permutation;
use super::transformation::UnimodularTransformation;

pub type Position = Vector3<f64>;
pub type AtomicSpecie = i32;

#[derive(Debug)]
pub struct Cell {
    pub lattice: Lattice,
    pub positions: Vec<Position>,
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

    pub fn num_atoms(&self) -> usize {
        self.positions.len()
    }

    pub fn transform(&self, transformation: &UnimodularTransformation) -> Self {
        let new_lattice = self.lattice.transform_unimodular(&transformation.linear);
        let pinv = transformation.linear_as_f64().try_inverse().unwrap();
        let new_positions = self
            .positions
            .iter()
            .map(|pos| pinv * (pos - transformation.origin_shift))
            .collect();
        Self {
            lattice: new_lattice,
            positions: new_positions,
            numbers: self.numbers.clone(),
        }
    }

    // Apply `trans`, which may increase the number of atoms in the cell.
    // Mapping from sites of the new cell to those of the original cell is also returned.
    // pub fn expand_transform(&self, transformation: &Transformation) -> (Self, SiteMapping) {
    //     unimplemented!()
    // }
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
        if identifier_mapping.contains_key(&uf.find(i)) {
            continue;
        }
        identifier_mapping.insert(uf.find(i), i);
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
    use crate::base::operation::Permutation;

    #[test]
    fn test_site_mapping_from_permutations() {
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
