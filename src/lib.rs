#[allow(unused_imports)]
#[macro_use]
extern crate approx;

pub mod base;
pub mod data;
pub mod identify;
pub mod math;
pub mod search;
pub mod symmetrize;

use nalgebra::Matrix3;

use crate::base::cell::Cell;
use crate::base::error::MoyoError;
use crate::base::operation::{AbstractOperations, Operations, Permutation};
use crate::base::tolerance::AngleTolerance;
use crate::base::transformation::Transformation;
use crate::data::hall_symbol_database::{HallNumber, Number};
use crate::data::setting::Setting;
use crate::identify::space_group::SpaceGroup;
use crate::search::primitive_cell::PrimitiveCellSearch;
use crate::search::symmetry_search::SymmetrySearch;

pub struct MoyoDataset {
    // Space-group type
    pub number: Number,
    pub hall_number: HallNumber,
    // Symmetry operations in the input cell
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
    // Site symmetry
    /// Spglib's `crystallographic_orbits` not `equivalent_atoms`
    pub orbits: Vec<usize>,
    // TODO: wyckoffs
    // TODO: site_symmetry_symbols
    // Standardized cell
    /// Transformation from the standardized cell to the input cell.
    pub std_transformation: Transformation,
    pub std_rotation: Matrix3<f64>,
    pub std_cell: Cell,
    // Standardized primitive cell
    /// Transformation from the standardized primitive cell to the input cell.
    pub std_prim_transformation: Transformation,
    pub std_prim_cell: Cell,
    pub mapping_to_std_prim: Vec<usize>,
}

impl MoyoDataset {
    pub fn new(
        cell: &Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        setting: Setting,
    ) -> Result<Self, MoyoError> {
        // Symmetry search
        let prim_cell_search = PrimitiveCellSearch::new(cell, symprec)?;
        let symmetry_search =
            SymmetrySearch::new(&prim_cell_search.primitive_cell, symprec, angle_tolerance)?;
        let prim_operations = AbstractOperations::from_operations(&symmetry_search.operations);

        // Space-group type identification
        let epsilon = symprec
            / prim_cell_search
                .primitive_cell
                .lattice
                .basis
                .determinant()
                .abs()
                .powf(1.0 / 3.0);
        let space_group = SpaceGroup::new(&prim_operations, setting, epsilon)?;

        // Ok(Self {
        //     number: space_group.number,
        //     hall_number: space_group.hall_number,
        // })
        unimplemented!();
    }

    pub fn num_operations(&self) -> usize {
        self.operations.num_operations()
    }
}
