extern crate nalgebra as na;
use na::base::Vector3;

use super::lattice::Lattice;
use super::transformation::UnimodularTransformation;

pub type Position = Vector3<f64>;
pub type AtomicSpecie = i32;
pub type SiteMapping = Vec<usize>;

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

#[cfg(test)]
mod tests {
    use std::panic;

    use nalgebra::{vector, Matrix3};

    use super::Cell;
    use crate::base::lattice::Lattice;

    #[test]
    fn test_mismatched_length() {
        let lattice = Lattice::new(Matrix3::<f64>::identity());
        let positions = vec![vector![0.0, 0.0, 0.0], vector![0.5, 0.5, 0.5]];
        let numbers = vec![1];

        let result = panic::catch_unwind(|| Cell::new(lattice, positions, numbers));
        assert!(result.is_err());
    }
}
