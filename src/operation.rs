use nalgebra::base::{Matrix3, Vector3};

use crate::cell::Cell;
use crate::lattice::Lattice;

pub type Rotation = Matrix3<i32>;
pub type Translation = Vector3<f64>;

pub type Permutation = Vec<usize>;

pub struct Operations {
    pub lattice: Lattice,
    //
    pub rotations: Vec<Rotation>,
    pub translations: Vec<Translation>,
}

pub struct CellWithOperations {
    pub cell: Cell,
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
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
}

impl CellWithOperations {
    pub fn new(cell: Cell, operations: Operations, permutations: Vec<Permutation>) -> Self {
        assert_relative_eq!(cell.lattice.basis, operations.lattice.basis);
        assert_eq!(operations.num_operations(), permutations.len());
        Self {
            cell,
            operations,
            permutations,
        }
    }

    pub fn num_operations(&self) -> usize {
        self.operations.num_operations()
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::{Operations, Translation};
    use crate::lattice::Lattice;

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
}
