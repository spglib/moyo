use std::collections::{HashSet, VecDeque};
use std::ops::Mul;

use nalgebra::base::{Matrix3, Vector3};

use super::lattice::Lattice;

/// Rotation matrix in a crystallographic basis
pub type Rotation = Matrix3<i32>;
/// Translation vector in a crystallographic basis
pub type Translation = Vector3<f64>;

#[derive(Debug, Clone)]
pub struct Operation {
    pub rotation: Rotation,
    pub translation: Translation,
}

impl Operation {
    pub fn new(rotation: Rotation, translation: Translation) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    /// Return rotation matrix in cartesian coordinates with respect to the given lattice
    pub fn cartesian_rotation(&self, lattice: &Lattice) -> Matrix3<f64> {
        let inv_basis = lattice.basis.try_inverse().unwrap();
        lattice.basis * self.rotation.map(|e| e as f64) * inv_basis
    }
}

impl Mul for Operation {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        // (r1, t1) * (r2, t2) = (r1 * r2, r1 * t2 + t1)
        let new_rotation = self.rotation * rhs.rotation;
        let new_translation = self.rotation.map(|e| e as f64) * rhs.translation + self.translation;
        Self::new(new_rotation, new_translation)
    }
}

pub type Rotations = Vec<Rotation>;
pub type Operations = Vec<Operation>;

pub fn project_rotations(operations: &Operations) -> Rotations {
    operations.iter().map(|ops| ops.rotation).collect()
}

#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct Permutation {
    pub mapping: Vec<usize>,
}

impl Permutation {
    pub fn new(mapping: Vec<usize>) -> Self {
        Self { mapping }
    }

    pub fn identity(size: usize) -> Self {
        Self::new((0..size).collect())
    }

    pub fn size(&self) -> usize {
        self.mapping.len()
    }

    pub fn apply(&self, i: usize) -> usize {
        self.mapping[i]
    }

    pub fn inverse(&self) -> Self {
        let mut inv = vec![0; self.size()];
        for (i, &j) in self.mapping.iter().enumerate() {
            inv[j] = i;
        }
        Self::new(inv)
    }
}

impl Mul for Permutation {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mapping = (0..self.size()).map(|i| self.apply(rhs.apply(i))).collect();
        Self::new(mapping)
    }
}

#[allow(dead_code)]
/// Used for testing
pub fn traverse(generators: &Vec<Rotation>) -> Vec<Rotation> {
    let mut queue = VecDeque::new();
    let mut visited = HashSet::new();
    let mut group = vec![];

    queue.push_back(Rotation::identity());

    while !queue.is_empty() {
        let element = queue.pop_front().unwrap();
        if visited.contains(&element) {
            continue;
        }
        visited.insert(element);
        group.push(element);

        for generator in generators {
            let product = element * generator;
            queue.push_back(product);
        }
    }

    group
}

#[cfg(test)]
mod tests {
    use std::vec;

    use nalgebra::matrix;

    use super::{Permutation, Translation};
    use crate::base::{lattice::Lattice, Operation};

    #[test]
    fn test_cartesian_rotations() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            -0.5, f64::sqrt(3.0) / 2.0, 0.0;
            0.0, 0.0, 1.0;
        ]);
        let operation = Operation::new(
            matrix![
                0, -1, 0;
                1, -1, 0;
                0, 0, 1;
            ],
            Translation::zeros(),
        );

        let actual = operation.cartesian_rotation(&lattice);
        let expect = matrix![
            -0.5, -f64::sqrt(3.0) / 2.0, 0.0;
            f64::sqrt(3.0) / 2.0, -0.5, 0.0;
            0.0, 0.0, 1.0;
        ];
        assert_relative_eq!(actual, expect);
    }

    #[test]
    fn test_permutation() {
        let permutation = Permutation::new(vec![1, 2, 0]);
        assert_eq!(permutation.apply(0), 1);
        assert_eq!(permutation.inverse(), Permutation::new(vec![2, 0, 1]));
        assert_eq!(
            permutation.clone() * permutation.inverse(),
            Permutation::identity(3)
        );
    }
}
