use nalgebra::vector;
use safer_ffi;
use safer_ffi::prelude::*;

use moyo::base::{AngleTolerance, Cell, Lattice};
use moyo::MoyoDataset;

use crate::base::MoyocOperation;
use crate::data::MoyocSetting;

#[derive_ReprC]
#[repr(C)]
#[derive(Debug)]
pub struct MoyocDataset {
    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    pub number: i32,
    pub hall_number: i32,
    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    pub operations: safer_ffi::Vec<MoyocOperation>,
    // // ------------------------------------------------------------------------
    // // Site symmetry
    // // ------------------------------------------------------------------------
    // pub orbits: Vec<usize>,
    // pub wyckoffs: Vec<char>,
    // pub site_symmetry_symbols: Vec<String>,
    // // ------------------------------------------------------------------------
    // // Standardized cell
    // // ------------------------------------------------------------------------
    // pub std_cell: Cell,
    // pub std_linear: Matrix3<f64>,
    // pub std_origin_shift: OriginShift,
    // pub std_rotation_matrix: Matrix3<f64>,
    // pub pearson_symbol: String,
    // // ------------------------------------------------------------------------
    // // Primitive standardized cell
    // // ------------------------------------------------------------------------
    // pub prim_std_cell: Cell,
    // pub prim_std_linear: Matrix3<f64>,
    // pub prim_std_origin_shift: OriginShift,
    // pub mapping_std_prim: Vec<usize>,
    // // ------------------------------------------------------------------------
    // // Final parameters
    // // ------------------------------------------------------------------------
    // pub symprec: f64,
    // pub angle_tolerance: AngleTolerance,
}

#[ffi_export]
pub fn moyoc_dataset(
    basis: *const [[f64; 3]; 3],
    positions: *const [f64; 3],
    numbers: *const i32,
    num_atoms: i32,
    symprec: f64,
    angle_tolerance: f64,
    setting: MoyocSetting,
) -> *const MoyocDataset {
    let basis = unsafe { &*basis };
    let lattice = Lattice::from_basis(*basis);
    let positions = unsafe {
        std::slice::from_raw_parts(positions, num_atoms as usize)
            .iter()
            .map(|&pos| vector![pos[0], pos[1], pos[2]])
            .collect::<Vec<_>>()
    };
    let numbers = unsafe {
        std::slice::from_raw_parts(numbers, num_atoms as usize)
            .iter()
            .map(|&number| number)
            .collect::<Vec<_>>()
    };
    let cell = Cell::new(lattice, positions, numbers);

    let angle_tolerance = if angle_tolerance < 0.0 {
        AngleTolerance::Default
    } else {
        AngleTolerance::Radian(angle_tolerance)
    };

    if let Ok(dataset) = MoyoDataset::new(&cell, symprec, angle_tolerance, setting.into()) {
        let dataset = MoyocDataset {
            // Identification
            number: dataset.number,
            hall_number: dataset.hall_number,
            // Symmetry operations in the input cell
            operations: dataset
                .operations
                .into_iter()
                .map(|ops| ops.into())
                .collect::<Vec<_>>()
                .into(),
        };
        return Box::into_raw(Box::new(dataset));
    } else {
        return std::ptr::null();
    }
}
