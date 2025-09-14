use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;
use moyo::MoyoDataset as Dataset;

use crate::base::MoyoCell;
use crate::data::MoyoSetting;

#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoDataset {
    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    pub number: i32,
    pub hall_number: i32,
    // pub hm_symbol: String,
    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    // pub operations: MoyoOperations,
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

impl From<Dataset> for MoyoDataset {
    fn from(dataset: Dataset) -> Self {
        Self {
            number: dataset.number,
            hall_number: dataset.hall_number,
        }
    }
}

#[no_mangle]
/// Create a dataset from the given cell.
/// hall_number: -1 means None
pub extern "C" fn moyo_dataset(
    basis: *const [[f64; 3]; 3],
    positions: *const [f64; 3],
    numbers: *const i32,
    num_atoms: i32,
    symprec: f64,
    angle_tolerance: f64,
    setting: MoyoSetting,
    hall_number: i32,
) -> *mut MoyoDataset {
    let cell = MoyoCell {
        basis: unsafe { *basis },
        positions,
        numbers,
        num_atoms,
    };
    let cell: Cell = (&cell).into();

    let angle_tolerance = if angle_tolerance < 0.0 {
        AngleTolerance::Default
    } else {
        AngleTolerance::Radian(angle_tolerance)
    };
    let setting = match setting {
        MoyoSetting::HallNumber => {
            if hall_number <= 0 {
                return std::ptr::null_mut();
            }
            Setting::HallNumber(hall_number)
        }
        MoyoSetting::Spglib => Setting::Spglib,
        MoyoSetting::Standard => Setting::Standard,
    };

    let dataset = Dataset::new(&cell, symprec, angle_tolerance, setting);
    let dataset_ptr = match dataset {
        Ok(dataset) => {
            let moyoc_dataset: MoyoDataset = dataset.into();
            Box::into_raw(Box::new(moyoc_dataset))
        }
        Err(_) => std::ptr::null_mut(),
    };
    dataset_ptr
}

#[no_mangle]
pub extern "C" fn free_moyo_dataset(dataset: *mut MoyoDataset) {}
