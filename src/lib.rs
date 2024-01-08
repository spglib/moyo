#[allow(unused_imports)]
#[macro_use]
extern crate approx;

pub mod base;
pub mod data;

mod identify;
mod math;
mod search;
mod symmetrize;

use crate::base::{
    orbits_from_permutations, AbstractOperations, AngleTolerance, Cell, MoyoError, OriginShift,
    Transformation,
};
use crate::data::{HallNumber, Number, Setting};
use crate::identify::SpaceGroup;
use crate::search::{PrimitiveCell, SymmetrySearch};

#[derive(Debug)]
pub struct MoyoDataset {
    // Space-group type
    pub number: Number,
    pub hall_number: HallNumber,
    // Symmetry operations in the input cell
    pub operations: AbstractOperations,
    /// Spglib's `crystallographic_orbits` not `equivalent_atoms`
    /// For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent and the last four atoms are equivalent.
    pub orbits: Vec<usize>,
    // Site symmetry
    // TODO: wyckoffs
    // TODO: site_symmetry_symbols
    // Standardized cell
    // Transformation from the standardized cell to the input cell.
    // TODO: pub std_transformation: Transformation,
    // TODO: pub std_rotation_matrix: Matrix3<f64>,
    // TODO: pub std_cell: Cell,
    // Standardized primitive cell
    // Transformation from the standardized primitive cell to the input cell.
    // TODO: pub std_prim_transformation: Transformation,
    // TODO: pub std_prim_cell: Cell,
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
        let prim_operations = AbstractOperations::from_operations(&symmetry_search.operations);

        // Symmetry in the input cell
        let operations = operations_in_cell(&prim_cell, &prim_operations);
        let orbits = orbits_in_cell(&prim_cell, &symmetry_search);

        // Space-group type identification
        let epsilon = symprec / prim_cell.cell.lattice.volume().powf(1.0 / 3.0);
        let space_group = SpaceGroup::new(&prim_operations, setting, epsilon)?;

        Ok(Self {
            number: space_group.number,
            hall_number: space_group.hall_number,
            operations,
            orbits,
        })
    }

    pub fn num_operations(&self) -> usize {
        self.operations.num_operations()
    }
}

fn operations_in_cell(
    prim_cell: &PrimitiveCell,
    prim_operations: &AbstractOperations,
) -> AbstractOperations {
    let mut rotations = vec![];
    let mut translations = vec![];
    let input_operations =
        prim_operations.transform(&Transformation::new(prim_cell.linear, OriginShift::zeros()));
    for t1 in prim_cell.translations.iter() {
        for (rotation, t2) in input_operations
            .rotations
            .iter()
            .zip(input_operations.translations.iter())
        {
            // (E, t1) (rotation, t2) = (rotation, t1 + t2)
            rotations.push(*rotation);
            let mut t12 = t1 + t2;
            t12 -= t12.map(|x| x.round());
            translations.push(t12);
        }
    }

    AbstractOperations::new(rotations, translations)
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
