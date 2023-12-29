use std::collections::{HashSet, VecDeque};
use std::ops::Mul;

use nalgebra::base::{Matrix3, Vector3};

use super::lattice::Lattice;

pub type Rotation = Matrix3<i32>;
pub type Translation = Vector3<f64>;

#[derive(Debug, PartialEq, Clone)]
pub struct Permutation {
    pub mapping: Vec<usize>,
}

#[derive(Debug)]
/// Symmetry operation without basis information
pub struct AbstractOperations {
    pub rotations: Vec<Rotation>,
    pub translations: Vec<Translation>,
}

impl AbstractOperations {
    pub fn new(rotations: Vec<Rotation>, translations: Vec<Translation>) -> Self {
        if translations.len() != rotations.len() {
            panic!("rotations and translations should be the same length");
        }
        Self {
            rotations,
            translations,
        }
    }

    pub fn num_operations(&self) -> usize {
        self.rotations.len()
    }

    /// (P, p)^-1 (W, w) (P, p)
    pub fn transform(&self, linear: &Matrix3<f64>, origin_shift: &Vector3<f64>) -> Self {
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in self.rotations.iter().zip(self.translations.iter()) {
            let (new_rotation, new_translation) =
                transform_operation(rotation, translation, linear, origin_shift);
            new_rotations.push(new_rotation);
            new_translations.push(new_translation);
        }
        Self::new(new_rotations, new_translations)
    }
}

fn transform_operation(
    rotation: &Rotation,
    translation: &Translation,
    linear: &Matrix3<f64>,
    origin_shift: &Vector3<f64>,
) -> (Rotation, Translation) {
    let linear_inv = linear.try_inverse().unwrap();
    let new_rotation = (linear_inv * rotation.map(|e| e as f64) * linear).map(|e| e.round() as i32);
    let new_translation =
        linear_inv * (rotation.map(|e| e as f64) * origin_shift + translation - origin_shift);
    (new_rotation, new_translation)
}

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

#[derive(Debug)]
/// Symmetry operation on a given basis
pub struct Operations {
    pub lattice: Lattice,
    //
    pub rotations: Vec<Rotation>,
    pub translations: Vec<Translation>,
}

impl Operations {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>, translations: Vec<Translation>) -> Self {
        if translations.len() != rotations.len() {
            panic!("rotations and translations should be the same length");
        }
        Self {
            lattice,
            rotations,
            translations,
        }
    }

    pub fn num_operations(&self) -> usize {
        self.rotations.len()
    }

    pub fn cartesian_rotations(&self) -> Vec<Matrix3<f64>> {
        let inv_basis = self.lattice.basis.try_inverse().unwrap();
        self.rotations
            .iter()
            .map(|r| self.lattice.basis * r.map(|e| e as f64) * inv_basis)
            // .map(|r| inv_basis * r.map(|e| e as f64) * self.lattice.basis)
            .collect()
    }

    pub fn transform(&self, linear: &Matrix3<f64>, origin_shift: &Vector3<f64>) -> Self {
        let new_lattice = self.lattice.transform(linear);
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in self.rotations.iter().zip(self.translations.iter()) {
            let (new_rotation, new_translation) =
                transform_operation(rotation, translation, linear, origin_shift);
            new_rotations.push(new_rotation);
            new_translations.push(new_translation);
        }
        Self::new(new_lattice, new_rotations, new_translations)
    }
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

#[cfg(test)]
mod tests {
    use std::vec;

    use nalgebra::matrix;

    use super::{Operations, Permutation, Translation};
    use crate::base::lattice::Lattice;

    #[test]
    fn test_cartesian_rotations() {
        let lattice = Lattice::new(matrix![
            1.0, -0.5, 0.0;
            0.0, f64::sqrt(3.0) / 2.0, 0.0;
            0.0, 0.0, 1.0;
        ]);
        let rotations = vec![matrix![
            0, -1, 0;
            1, -1, 0;
            0, 0, 1;
        ]];
        let translations = vec![Translation::zeros()];

        let operations = Operations::new(lattice, rotations, translations);

        let actual = operations.cartesian_rotations()[0];
        let expect = matrix![
            -0.5, -f64::sqrt(3.0) / 2.0, 0.0;
            f64::sqrt(3.0) / 2.0, -0.5, 0.0;
            0.0, 0.0, 1.0;
        ];
        assert_relative_eq!(actual, expect);
        assert_eq!(operations.num_operations(), 1)
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
