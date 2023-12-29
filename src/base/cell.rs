extern crate nalgebra as na;
use na::base::Vector3;

use super::lattice::Lattice;
use super::transformation::Transformation;

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

    pub fn transform(&self, transformation: &Transformation) -> Self {
        if transformation.size != 1 {
            panic!("transformation matrix should have determinant of -1 or 1");
        }
        let new_lattice = self
            .lattice
            .transform(&transformation.trans_mat.map(|e| e as f64));
        let pinv = transformation.trans_mat_as_f64().try_inverse().unwrap();
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
