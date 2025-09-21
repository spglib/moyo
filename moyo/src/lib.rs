/*!
# moyo

**moyo** is a fast and robust crystal symmetry finder.

## Using **moyo**

Simply add the following to your `Cargo.toml` file:

```ignore
[dependencies]
// TODO: replace the * with the latest version
moyo = "*"
```

## Examples

The basic usage of **moyo** is to create a [`moyo::base::Cell`](Cell) representing a crystal structure, and then create a [`moyo::MoyoDataset`](MoyoDataset) from the [`moyo::base::Cell`](Cell).
The [`moyo::MoyoDataset`](MoyoDataset) contains symmetry information of the input crystal structure: for example, the space group number, symmetry operations, and standardized cell.

Magnetic symmetry is also supported in **moyo**.
Magnetic moments are represented by a struct implementing the [`moyo::base::MagneticMoment`](MagneticMoment) trait: for example, [`moyo::base::Collinear`](Collinear) or [`moyo::base::NonCollinear`](NonCollinear).
Magnetic cell is represented by a [`moyo::base::MagneticCell`](MagneticCell) struct.
The [`moyo::MoyoMagneticDataset`](MoyoMagneticDataset) contains magnetic symmetry information of the input magnetic cell: for example, the magnetic space-group type, magnetic symmetry operations, and standardized magnetic cell.

```
use nalgebra::{matrix, vector, Matrix3, Vector3};
use moyo::{MoyoDataset, MoyoMagneticDataset};
use moyo::base::{Cell, MagneticCell, AngleTolerance, Lattice, NonCollinear, RotationMagneticMomentAction};
use moyo::data::Setting;

let lattice = Lattice::new(matrix![
    4.603, 0.0, 0.0;
    0.0, 4.603, 0.0;
    0.0, 0.0, 2.969;
]);
let x_4f = 0.3046;
let positions = vec![
    Vector3::new(0.0, 0.0, 0.0),                // Ti(2a)
    Vector3::new(0.5, 0.5, 0.5),                // Ti(2a)
    Vector3::new(x_4f, x_4f, 0.0),              // O(4f)
    Vector3::new(-x_4f, -x_4f, 0.0),            // O(4f)
    Vector3::new(-x_4f + 0.5, x_4f + 0.5, 0.5), // O(4f)
    Vector3::new(x_4f + 0.5, -x_4f + 0.5, 0.5), // O(4f)
];
let numbers = vec![0, 0, 1, 1, 1, 1];
let cell = Cell::new(lattice.clone(), positions.clone(), numbers.clone());

let symprec = 1e-4;
let angle_tolerance = AngleTolerance::Default;
let setting = Setting::Standard;

let dataset = MoyoDataset::new(&cell, symprec, angle_tolerance, setting).unwrap();
assert_eq!(dataset.number, 136);  // P4_2/mnm

let magnetic_moments = vec![
    NonCollinear(vector![0.0, 0.0, 0.7]),
    NonCollinear(vector![0.0, 0.0, -0.7]),
    NonCollinear(vector![0.0, 0.0, 0.0]),
    NonCollinear(vector![0.0, 0.0, 0.0]),
    NonCollinear(vector![0.0, 0.0, 0.0]),
    NonCollinear(vector![0.0, 0.0, 0.0]),
];
let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magnetic_moments);

let action = RotationMagneticMomentAction::Axial;

let magnetic_dataset = MoyoMagneticDataset::new(&magnetic_cell, symprec, angle_tolerance, None, action).unwrap();
assert_eq!(magnetic_dataset.uni_number, 1159);  // BNS 136.499
```

## Features
- Support most of the symmetry search functionality in Spglib
- Primitive cell search
- Symmetry operation search
- Space-group type identification
- Wyckoff position assignment
- Crystal structure symmetrization
- Magnetic space group support

*/
#[allow(unused_imports)]
#[macro_use]
extern crate approx;

pub mod base;
pub mod data;
pub mod identify;
pub mod math;
pub mod search; // Public for benchmarking
pub mod utils;

mod symmetrize;

use crate::base::{
    AngleTolerance, Cell, MagneticCell, MagneticMoment, MagneticOperations, MoyoError, Operations,
    OriginShift, RotationMagneticMomentAction,
};
use crate::data::{
    HallNumber, Number, Setting, UNINumber, arithmetic_crystal_class_entry, hall_symbol_entry,
};
use crate::identify::{MagneticSpaceGroup, SpaceGroup};
use crate::search::{
    iterative_magnetic_symmetry_search, iterative_symmetry_search,
    magnetic_operations_in_magnetic_cell, operations_in_cell,
};
use crate::symmetrize::{StandardizedCell, StandardizedMagneticCell, orbits_in_cell};

use nalgebra::Matrix3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, Serialize)]
/// A dataset containing symmetry information of the input crystal structure.
pub struct MoyoDataset {
    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    /// Space group number.
    pub number: Number,
    /// Hall symbol number.
    pub hall_number: HallNumber,
    /// Hermann-Mauguin symbol in short notation (e.g., "Fd-3m" for space group 227).
    pub hm_symbol: String,
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
    /// Pearson symbol for standardized cell
    pub pearson_symbol: String,
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
    /// Create a new [`MoyoDataset`] from the input cell, `cell`.
    /// `symprec` and `angle_tolerance` control the tolerances for searching symmetry operations.
    /// `setting` determines the preference for the "standardized" setting of a detected space-group type.
    /// If the search fails, [`MoyoError`] is returned.
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
        let std_cell = StandardizedCell::new(
            &prim_cell.cell,
            &symmetry_search.operations,
            &symmetry_search.permutations,
            &space_group,
            symprec,
            epsilon,
        )?;

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

        // Pearson symbol
        let hall_symbol = hall_symbol_entry(space_group.hall_number).unwrap();
        let arithmetic_entry =
            arithmetic_crystal_class_entry(hall_symbol.arithmetic_number).unwrap();
        let bravais_class = arithmetic_entry.bravais_class;
        let pearson_symbol = format!("{}{}", bravais_class.to_string(), std_cell.cell.num_atoms());

        Ok(Self {
            // Space-group type
            number: space_group.number,
            hall_number: space_group.hall_number,
            hm_symbol: hall_symbol.hm_short.to_string(),
            // Symmetry operations in the input cell
            operations,
            // Standardized cell
            std_cell: std_cell.cell,
            std_linear,
            std_origin_shift,
            std_rotation_matrix: std_cell.rotation_matrix,
            pearson_symbol,
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

    /// Return the number of symmetry operations in the input cell.
    pub fn num_operations(&self) -> usize {
        self.operations.len()
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MoyoMagneticDataset<M: MagneticMoment> {
    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------
    pub uni_number: UNINumber,
    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Magnetic symmetry operations in the input cell.
    pub magnetic_operations: MagneticOperations,
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom in the **input** magnetic cell.
    /// For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent and the last four atoms are equivalent to each other.
    pub orbits: Vec<usize>,
    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Standardized magnetic cell
    pub std_mag_cell: MagneticCell<M>,
    /// Linear part of transformation from the input magnetic cell to the standardized one.
    pub std_linear: Matrix3<f64>,
    /// Origin shift of transformation from the input magnetic cell to the standardized one.
    pub std_origin_shift: OriginShift,
    /// Rigid rotation
    pub std_rotation_matrix: Matrix3<f64>,
    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    pub prim_std_mag_cell: MagneticCell<M>,
    /// Linear part of transformation from the input magnetic cell to the primitive standardized magnetic cell.
    pub prim_std_linear: Matrix3<f64>,
    /// Origin shift of transformation from the input magnetic cell to the primitive standardized magnetic cell.
    pub prim_std_origin_shift: OriginShift,
    /// Mapping sites in the input magnetic cell to those in the primitive standardized magnetic cell.
    /// The `i`th atom in the input magnetic cell is mapped to the `mapping_to_std_prim[i]`th atom in the primitive standardized magnetic cell.
    pub mapping_std_prim: Vec<usize>,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    pub angle_tolerance: AngleTolerance,
    /// Actually used `mag_symprec` in iterative symmetry search.
    pub mag_symprec: f64,
}

impl<M: MagneticMoment> MoyoMagneticDataset<M> {
    pub fn new(
        magnetic_cell: &MagneticCell<M>,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        mag_symprec: Option<f64>,
        action: RotationMagneticMomentAction,
    ) -> Result<Self, MoyoError> {
        let (prim_mag_cell, magnetic_symmetry_search, symprec, angle_tolerance, mag_symprec) =
            iterative_magnetic_symmetry_search(
                magnetic_cell,
                symprec,
                angle_tolerance,
                mag_symprec,
                action,
            )?;
        let magnetic_operations = magnetic_operations_in_magnetic_cell(
            &prim_mag_cell,
            &magnetic_symmetry_search.magnetic_operations,
        );

        // Magnetic space-group type identification
        let epsilon = symprec
            / prim_mag_cell
                .magnetic_cell
                .cell
                .lattice
                .volume()
                .powf(1.0 / 3.0);
        let magnetic_space_group =
            MagneticSpaceGroup::new(&magnetic_symmetry_search.magnetic_operations, epsilon)?;

        // Standardized magnetic cell
        let std_mag_cell = StandardizedMagneticCell::new(
            &prim_mag_cell,
            &magnetic_symmetry_search,
            &magnetic_space_group,
            symprec,
            mag_symprec,
            epsilon,
            action,
        )?;

        // Site symmetry
        // StandardizedMagneticCell.prim_mag_cell and prim_mag_cell have the same site order
        let mapping_std_prim = prim_mag_cell.site_mapping.clone();
        let orbits = orbits_in_cell(
            prim_mag_cell.magnetic_cell.num_atoms(),
            &magnetic_symmetry_search.permutations,
            &mapping_std_prim,
        );

        // magnetic_cell <-(prim_mag_cell.linear, 0)- prim_mag_cell.magnetic_cell -(std_mag_cell.transformation)-> std_mag_cell.mag_cell
        // (std_linear, std_origin_shift) = (prim_mag_cell.linear^-1, 0) * std_mag_cell.transformation
        let prim_mag_cell_linear_inv = prim_mag_cell
            .linear
            .map(|e| e as f64)
            .try_inverse()
            .unwrap();
        let std_linear = prim_mag_cell_linear_inv * std_mag_cell.transformation.linear_as_f64();
        let std_origin_shift = prim_mag_cell_linear_inv * std_mag_cell.transformation.origin_shift;

        // (prim_std_linear, prim_std_origin_shift) = (prim_mag_cell.linear^-1, 0) * std_mag_cell.prim_transformation
        let prim_std_linear =
            prim_mag_cell_linear_inv * std_mag_cell.prim_transformation.linear_as_f64();
        let prim_std_origin_shift =
            prim_mag_cell_linear_inv * std_mag_cell.prim_transformation.origin_shift;

        Ok(Self {
            // Magnetic space-group type
            uni_number: magnetic_space_group.uni_number,
            // Magnetic symmetry operations in the input cell
            magnetic_operations,
            // Site symmetry
            orbits,
            // Standardized magnetic cell
            std_mag_cell: std_mag_cell.mag_cell,
            std_linear,
            std_origin_shift,
            std_rotation_matrix: std_mag_cell.rotation_matrix,
            // Primitive standardized magnetic cell
            prim_std_mag_cell: std_mag_cell.prim_mag_cell,
            prim_std_linear,
            prim_std_origin_shift,
            mapping_std_prim,
            // Final parameters
            symprec,
            angle_tolerance,
            mag_symprec,
        })
    }

    /// Return the number of magnetic symmetry operations in the input magnetic cell.
    pub fn num_magnetic_operations(&self) -> usize {
        self.magnetic_operations.len()
    }
}
