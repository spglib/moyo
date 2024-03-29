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

use crate::base::{AngleTolerance, Cell, MoyoError, Operations, OriginShift, Transformation};
use crate::data::{HallNumber, Number, Setting};
use crate::identify::SpaceGroup;
use crate::search::{PrimitiveCell, PrimitiveSymmetrySearch};
use crate::symmetrize::{orbits_in_cell, StandardizedCell};

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
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// Spglib's `crystallographic_orbits` not `equivalent_atoms`
    /// The `i`th atom in the input cell is equivalent to the `orbits[i]`th atom in the **input** cell.
    /// For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent and the last four atoms are equivalent to each other.
    pub orbits: Vec<usize>,
    /// Wyckoff letters for each site in the input cell.
    pub wyckoffs: Vec<char>,
    /// Site symmetry symbols for each site in the input cell.
    /// The orientation of the site symmetry is w.r.t. the standardized cell.
    pub site_symmetry_symbols: Vec<String>,
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
    /// Mapping sites in the input cell to those in the primitive standardized cell.
    /// The `i`th atom in the input cell is mapped to the `mapping_to_std_prim[i]`th atom in the primitive standardized cell.
    pub mapping_std_prim: Vec<usize>,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    // pub symprec: f64,
    // pub angle_tolerance: AngleTolerance,
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
        let symmetry_search =
            PrimitiveSymmetrySearch::new(&prim_cell.cell, symprec, angle_tolerance)?;
        let operations = operations_in_cell(&prim_cell, &symmetry_search.operations);

        // Space-group type identification
        let epsilon = symprec / prim_cell.cell.lattice.volume().powf(1.0 / 3.0);
        let space_group = SpaceGroup::new(&symmetry_search.operations, setting, epsilon)?;

        // Standardized cell
        let std_cell = StandardizedCell::new(&prim_cell, &symmetry_search, &space_group, symprec)?;

        // site symmetry
        let orbits = orbits_in_cell(
            prim_cell.cell.num_atoms(),
            &symmetry_search.permutations,
            &prim_cell.site_mapping,
        );
        // StandardizedCell.prim_cell and prim_cell have the same site order
        let mapping_std_prim = prim_cell.site_mapping.clone();
        let mut std_prim_wyckoffs = vec![None; prim_cell.cell.num_atoms()];
        for (i, wyckoff) in std_cell.wyckoffs.iter().enumerate() {
            let j = std_cell.site_mapping[i];
            if std_prim_wyckoffs[j].is_none() {
                std_prim_wyckoffs[j] = Some(wyckoff.clone());
            }
        }
        let wyckoffs: Option<Vec<_>> = mapping_std_prim
            .iter()
            .map(|&i| std_prim_wyckoffs[i].clone())
            .collect();
        let wyckoffs = wyckoffs.ok_or(MoyoError::WyckoffPositionAssignmentError)?;

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
            // Standardized cell
            std_cell: std_cell.cell,
            std_linear,
            std_origin_shift,
            std_rotation_matrix: std_cell.rotation_matrix,
            // Primitive standardized cell
            prim_std_cell: std_cell.prim_cell,
            prim_std_linear,
            prim_std_origin_shift,
            mapping_std_prim,
            // Site symmetry
            orbits,
            wyckoffs: wyckoffs.iter().map(|w| w.letter).collect(),
            site_symmetry_symbols: wyckoffs
                .iter()
                .map(|w| w.site_symmetry.to_string())
                .collect(),
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
