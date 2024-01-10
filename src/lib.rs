#[allow(unused_imports)]
#[macro_use]
extern crate approx;

pub mod base;
pub mod data;

mod identify;
mod math;
mod search;
mod symmetrize;

use nalgebra::Matrix3;

use crate::base::{
    orbits_from_permutations, AngleTolerance, Cell, MoyoError, Operations, OriginShift,
    Transformation,
};
use crate::data::{HallNumber, Number, Setting};
use crate::identify::SpaceGroup;
use crate::search::{PrimitiveCell, SymmetrySearch};
use crate::symmetrize::StandardizedCell;

#[derive(Debug)]
pub struct MoyoDataset {
    // ------------------------------------------------------------------------
    // Space-group type
    // ------------------------------------------------------------------------
    pub number: Number,
    pub hall_number: HallNumber,
    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    pub operations: Operations,
    /// Spglib's `crystallographic_orbits` not `equivalent_atoms`
    /// For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent and the last four atoms are equivalent.
    pub orbits: Vec<usize>,
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    // TODO: wyckoffs
    // TODO: site_symmetry_symbols
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    pub std_cell: Cell,
    /// Linear part of transformation from the input cell to the standardized cell.
    pub std_linear: Matrix3<f64>,
    /// Origin shift of transformation from the input cell to the standardized cell.
    pub std_origin_shift: OriginShift,
    /// Rigid rotation
    pub std_rotation_matrix: Matrix3<f64>,
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    pub prim_std_cell: Cell,
    /// Linear part of transformation from the input cell to the primitive standardized cell.
    pub prim_std_linear: Matrix3<f64>,
    /// Origin shift of transformation from the input cell to the primitive standardized cell.
    pub prim_std_origin_shift: OriginShift,
    // TODO: pub mapping_to_std_prim: Vec<usize>,
    // TODO: pub std_prim_permutations: Vec<Permutation>,
}

impl MoyoDataset {
    pub fn new(
        cell: &Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        setting: Setting,
    ) -> Result<Self, MoyoError> {
        // Symmetry search
        let prim_cell = PrimitiveCell::new(cell, symprec)?;
        let symmetry_search = SymmetrySearch::new(&prim_cell.cell, symprec, angle_tolerance)?;

        // Symmetry in the input cell
        let operations = operations_in_cell(&prim_cell, &symmetry_search.operations);
        let orbits = orbits_in_cell(&prim_cell, &symmetry_search);

        // Space-group type identification
        let epsilon = symprec / prim_cell.cell.lattice.volume().powf(1.0 / 3.0);
        let space_group = SpaceGroup::new(&symmetry_search.operations, setting, epsilon)?;

        // Standardized cell
        let std_cell = StandardizedCell::new(&prim_cell, &space_group)?;

        // cell <-(prim_cell.linear, 0)- prim_cell.cell -(std_cell.transformation)-> std_cell.cell
        // (std_linear, std_origin_shift) = (prim_cell.linear^-1, 0) * std_cell.transformation
        let prim_cell_linear_inv = prim_cell.linear.map(|e| e as f64).try_inverse().unwrap();
        let std_linear = prim_cell_linear_inv * std_cell.transformation.linear_as_f64();
        let std_origin_shift = prim_cell_linear_inv * std_cell.transformation.origin_shift;

        // (prim_std_linear, prim_std_origin_shift) = (prim_cell.linear^-1, 0) * std_cell.prim_transformation
        let prim_std_linear = prim_cell_linear_inv * std_cell.prim_transformation.linear_as_f64();
        let prim_std_origin_shift =
            prim_cell_linear_inv * std_cell.prim_transformation.origin_shift;

        Ok(Self {
            // Space-group type
            number: space_group.number,
            hall_number: space_group.hall_number,
            // Symmetry operations in the input cell
            operations,
            orbits,
            // Standardized cell
            std_cell: std_cell.cell,
            std_linear,
            std_origin_shift,
            std_rotation_matrix: std_cell.rotation_matrix,
            // Primitive standardized cell
            prim_std_cell: std_cell.prim_cell,
            prim_std_linear,
            prim_std_origin_shift,
        })
    }

    pub fn num_operations(&self) -> usize {
        self.operations.num_operations()
    }
}

fn operations_in_cell(prim_cell: &PrimitiveCell, prim_operations: &Operations) -> Operations {
    let mut rotations = vec![];
    let mut translations = vec![];
    let input_operations =
        Transformation::from_linear(prim_cell.linear).transform_operations(prim_operations);
    for t1 in prim_cell.translations.iter() {
        for (rotation, t2) in input_operations
            .rotations
            .iter()
            .zip(input_operations.translations.iter())
        {
            // (E, t1) (rotation, t2) = (rotation, t1 + t2)
            rotations.push(*rotation);
            let t12 = (t1 + t2).map(|e| e % 1.);
            translations.push(t12);
        }
    }

    Operations::new(rotations, translations)
}

fn orbits_in_cell(prim_cell: &PrimitiveCell, symmetry_search: &SymmetrySearch) -> Vec<usize> {
    // prim_site_mapping: [prim_num_atoms] -> [prim_num_atoms]
    let prim_orbits =
        orbits_from_permutations(prim_cell.cell.num_atoms(), &symmetry_search.permutations);

    let num_atoms = prim_cell.site_mapping.len();
    (0..num_atoms)
        .map(|i| prim_orbits[prim_cell.site_mapping[i]])
        .collect::<Vec<_>>()
}
