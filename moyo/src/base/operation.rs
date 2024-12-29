use std::collections::{HashSet, VecDeque};
use std::ops::Mul;

use nalgebra::base::{Matrix3, Vector3};

use super::lattice::Lattice;

/// Rotation matrix in a crystallographic basis
pub type Rotation = Matrix3<i32>;
pub type CartesianRotation = Matrix3<f64>;
/// Translation vector in a crystallographic basis
pub type Translation = Vector3<f64>;
/// Time reversal operation
pub type TimeReversal = bool;

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
    pub fn cartesian_rotation(&self, lattice: &Lattice) -> CartesianRotation {
        let inv_basis = lattice.basis.try_inverse().unwrap();
        lattice.basis * self.rotation.map(|e| e as f64) * inv_basis
    }

    pub fn identity() -> Self {
        Self::new(Rotation::identity(), Translation::zeros())
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

#[derive(Debug, Clone)]
pub struct MagneticOperation {
    pub operation: Operation,
    pub time_reversal: TimeReversal,
}

impl MagneticOperation {
    pub fn new(rotation: Rotation, translation: Translation, time_reversal: TimeReversal) -> Self {
        let operation = Operation::new(rotation, translation);
        Self::from_operation(operation, time_reversal)
    }

    pub fn from_operation(operation: Operation, time_reversal: TimeReversal) -> Self {
        Self {
            operation,
            time_reversal,
        }
    }

    pub fn identity() -> Self {
        Self::from_operation(Operation::identity(), false)
    }
}

pub type Rotations = Vec<Rotation>;
pub type Operations = Vec<Operation>;
pub type MagneticOperations = Vec<MagneticOperation>;

pub fn project_rotations(operations: &Operations) -> Rotations {
    operations.iter().map(|ops| ops.rotation).collect()
}

pub fn project_operations(magnetic_operations: &MagneticOperations) -> Operations {
    magnetic_operations
        .iter()
        .map(|mops| mops.operation.clone())
        .collect()
}

#[allow(dead_code)]
/// Used for testing
pub fn traverse(generators: &Rotations) -> Rotations {
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
    use nalgebra::matrix;

    use super::Translation;
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
}
