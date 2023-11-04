extern crate nalgebra as na;
use na::base::Vector3;

use crate::lattice::Lattice;
use crate::transformation::Transformation;

pub type Position = Vector3<f64>;
pub type AtomicSpecie = i32;

#[derive(Debug)]
pub struct Cell {
    pub lattice: Lattice,
    pub positions: Vec<Position>,
    pub numbers: Vec<AtomicSpecie>,
    //
    pub num_atoms: usize,
}

impl Cell {
    pub fn new(lattice: Lattice, positions: Vec<Position>, numbers: Vec<AtomicSpecie>) -> Self {
        let num_atoms = positions.len();
        if numbers.len() != num_atoms {
            panic!("positions and numbers should be the same length");
        }
        Self {
            lattice,
            positions,
            numbers,
            num_atoms,
        }
    }

    pub fn transform(&self, trans: &Transformation) -> Self {
        if trans.size == 1 {
            self.transform_unimodular(trans)
        } else {
            unimplemented!("transformation matrix should have determinant of -1 or 1");
        }
    }

    fn transform_unimodular(&self, trans: &Transformation) -> Self {
        let new_lattice = self.lattice.transform(&trans);
        let pinv = trans.trans_mat_as_f64().try_inverse().unwrap();
        let new_positions = self
            .positions
            .iter()
            .map(|pos| pinv * (pos - trans.origin_shift))
            .collect();
        Self {
            lattice: new_lattice,
            positions: new_positions,
            numbers: self.numbers.clone(),
            num_atoms: self.num_atoms,
        }
    }
}
