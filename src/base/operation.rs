use std::collections::{HashSet, VecDeque};
use std::ops::Mul;

use nalgebra::base::{Matrix3, Vector3};

use super::lattice::Lattice;
use super::transformation::{
    Linear, OriginShift, Transformation, UnimodularLinear, UnimodularTransformation,
};

pub type Rotation = Matrix3<i32>;
pub type Translation = Vector3<f64>;

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

    pub fn from_operations(operations: &Operations) -> Self {
        Self::new(
            operations.operations.rotations.clone(),
            operations.operations.translations.clone(),
        )
    }

    pub fn num_operations(&self) -> usize {
        self.rotations.len()
    }

    /// (P, p)^-1 (W, w) (P, p)
    /// This function may decrease the number of operations if the transformation is not compatible with an operation.
    pub fn transform(&self, transformation: &Transformation) -> Self {
        let linear_inv = transformation.linear.try_inverse().unwrap();
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in self.rotations.iter().zip(self.translations.iter()) {
            if let Some((new_rotation, new_translation)) = transform_operation(
                rotation,
                translation,
                &transformation.linear,
                &linear_inv,
                &transformation.origin_shift,
            ) {
                new_rotations.push(new_rotation);
                new_translations.push(new_translation);
            }
        }
        Self::new(new_rotations, new_translations)
    }

    pub fn transform_unimodular(&self, transformation: &UnimodularTransformation) -> Self {
        let linear_inv = transformation
            .linear_as_f64()
            .try_inverse()
            .unwrap()
            .map(|e| e.round() as i32);
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in self.rotations.iter().zip(self.translations.iter()) {
            let (new_rotation, new_translation) = transform_operation_unimodular(
                rotation,
                translation,
                &transformation.linear,
                &linear_inv,
                &transformation.origin_shift,
            );
            new_rotations.push(new_rotation);
            new_translations.push(new_translation);
        }
        Self::new(new_rotations, new_translations)
    }
}

fn transform_operation(
    rotation: &Rotation,
    translation: &Translation,
    linear: &Linear,
    linear_inv: &Linear,
    origin_shift: &OriginShift,
) -> Option<(Rotation, Translation)> {
    let new_rotation = (linear_inv * rotation.map(|e| e as f64) * linear).map(|e| e.round() as i32);

    // Check if `new_rotation` is an integer matrix
    let recovered =
        (linear * new_rotation.map(|e| e as f64) * linear_inv).map(|e| e.round() as i32);
    if recovered != *rotation {
        return None;
    }

    let new_translation =
        linear_inv * (rotation.map(|e| e as f64) * origin_shift + translation - origin_shift);
    Some((new_rotation, new_translation))
}

fn transform_operation_unimodular(
    rotation: &Rotation,
    translation: &Translation,
    linear: &UnimodularLinear,
    linear_inv: &UnimodularLinear,
    origin_shift: &OriginShift,
) -> (Rotation, Translation) {
    let new_rotation = linear_inv * rotation * linear;
    let new_translation = linear_inv.map(|e| e as f64)
        * (rotation.map(|e| e as f64) * origin_shift + translation - origin_shift);
    (new_rotation, new_translation)
}

#[derive(Debug)]
/// Symmetry operation on a given basis
pub struct Operations {
    pub lattice: Lattice,
    pub operations: AbstractOperations,
}

impl Operations {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>, translations: Vec<Translation>) -> Self {
        Self {
            lattice,
            operations: AbstractOperations::new(rotations, translations),
        }
    }

    pub fn num_operations(&self) -> usize {
        self.operations.rotations.len()
    }

    pub fn cartesian_rotations(&self) -> Vec<Matrix3<f64>> {
        let inv_basis = self.lattice.basis.try_inverse().unwrap();
        self.operations
            .rotations
            .iter()
            .map(|r| self.lattice.basis * r.map(|e| e as f64) * inv_basis)
            .collect()
    }

    pub fn transform(&self, transformation: &Transformation) -> Self {
        let new_lattice = self.lattice.transform(&transformation.linear);
        let new_operations = self.operations.transform(transformation);
        Self::new(
            new_lattice,
            new_operations.rotations,
            new_operations.translations,
        )
    }

    pub fn transform_unimodular(&self, transformation: &UnimodularTransformation) -> Self {
        let new_lattice = self.lattice.transform_unimodular(&transformation.linear);
        let new_operations = self.operations.transform_unimodular(transformation);
        Self::new(
            new_lattice,
            new_operations.rotations,
            new_operations.translations,
        )
    }
}

#[derive(Debug, PartialEq, Clone)]
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

    use nalgebra::{matrix, Matrix3};

    use super::{Operations, Permutation, Translation};
    use crate::base::lattice::Lattice;
    use crate::base::transformation::Transformation;

    #[test]
    fn test_incompatible_transformation() {
        let transformation = Transformation::from_linear(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 2.0;
        ]);
        // threefold rotation
        let operations = Operations::new(
            Lattice::new(Matrix3::<f64>::identity()),
            vec![matrix![
                0, 0, 1;
                1, 0, 0;
                0, 1, 0;
            ]],
            vec![Translation::zeros()],
        );
        assert_eq!(operations.transform(&transformation).num_operations(), 0);
    }

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
