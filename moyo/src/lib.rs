#[allow(unused_imports)]
#[macro_use]
extern crate approx;

pub mod base;
pub mod data;
pub mod identify;
pub mod math;
pub mod search;
pub mod symmetrize;

pub use base::{
    AngleTolerance, Cell, Lattice, MoyoError, Operations, OriginShift, Rotation, Translation,
};
pub use data::{HallNumber, Number, Setting};

use nalgebra::Matrix3;

use crate::base::{ToleranceHandler, Transformation};
use crate::identify::SpaceGroup;
use crate::search::{PrimitiveCell, PrimitiveSymmetrySearch};
use crate::symmetrize::{orbits_in_cell, StandardizedCell};

#[derive(Debug)]
pub struct MoyoDataset {
    // ------------------------------------------------------------------------
    // Space-group type
    // ------------------------------------------------------------------------
    /// Space group number.
    pub number: Number,
    /// Hall symbol number.
    pub hall_number: HallNumber,
    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Symmetry operations in the input cell.
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
    /// Standardized cell
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
    /// Primitive standardized cell
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
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    pub angle_tolerance: AngleTolerance,
}

impl MoyoDataset {
    /// Create a new `MoyoDataset` from the input cell.
    ///
    /// # Examples
    /// ```
    /// use nalgebra::{matrix, vector, Matrix3, Vector3};
    /// use moyo::{MoyoDataset, Cell, AngleTolerance, Setting, Lattice};
    /// let lattice = Lattice::new(matrix![
    ///     4.603, 0.0, 0.0;
    ///     0.0, 4.603, 0.0;
    ///     0.0, 0.0, 2.969;
    /// ]);
    /// let x_4f = 0.3046;
    /// let positions = vec![
    ///     Vector3::new(0.0, 0.0, 0.0),                // Ti(2a)
    ///     Vector3::new(0.5, 0.5, 0.5),                // Ti(2a)
    ///     Vector3::new(x_4f, x_4f, 0.0),              // O(4f)
    ///     Vector3::new(-x_4f, -x_4f, 0.0),            // O(4f)
    ///     Vector3::new(-x_4f + 0.5, x_4f + 0.5, 0.5), // O(4f)
    ///     Vector3::new(x_4f + 0.5, -x_4f + 0.5, 0.5), // O(4f)
    /// ];
    /// let numbers = vec![0, 0, 1, 1, 1, 1];
    /// let cell = Cell::new(lattice, positions, numbers);
    /// let symprec = 1e-5;
    /// let angle_tolerance = AngleTolerance::Default;
    /// let setting = Setting::Standard;
    /// let dataset = MoyoDataset::new(&cell, symprec, angle_tolerance, setting).unwrap();
    /// assert_eq!(dataset.number, 136);  // P4_2/mnm
    /// ```
    pub fn new(
        cell: &Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        setting: Setting,
    ) -> Result<Self, MoyoError> {
        let (prim_cell, symmetry_search, symprec, angle_tolerance) =
            iterative_symmetry_search(cell, symprec, angle_tolerance)?;
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
            // Final parameters
            symprec,
            angle_tolerance,
        })
    }

    pub fn num_operations(&self) -> usize {
        self.operations.num_operations()
    }
}

const MAX_SYMMETRY_SEARCH_TRIALS: usize = 16;

fn iterative_symmetry_search(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<(PrimitiveCell, PrimitiveSymmetrySearch, f64, AngleTolerance), MoyoError> {
    let mut tolerance_handler = ToleranceHandler::new(symprec, angle_tolerance);

    for _ in 0..MAX_SYMMETRY_SEARCH_TRIALS {
        match PrimitiveCell::new(cell, tolerance_handler.symprec) {
            Ok(prim_cell) => {
                match PrimitiveSymmetrySearch::new(
                    &prim_cell.cell,
                    tolerance_handler.symprec,
                    tolerance_handler.angle_tolerance,
                ) {
                    Ok(symmetry_search) => {
                        return Ok((
                            prim_cell,
                            symmetry_search,
                            tolerance_handler.symprec,
                            tolerance_handler.angle_tolerance,
                        ));
                    }
                    Err(err) => tolerance_handler.update(err),
                }
            }
            Err(err) => tolerance_handler.update(err),
        }
    }
    Err(MoyoError::PrimitiveSymmetrySearchError)
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
