use core::ffi::c_char;
use std::ffi::CString;

use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;
use moyo::MoyoDataset as Dataset;

use crate::base::{free_moyo_operations, MoyoCell, MoyoOperations};
use crate::data::MoyoSetting;

#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoDataset {
    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    /// Space group number.
    pub number: i32,
    /// Hall symbol number.
    pub hall_number: i32,
    /// Hermann-Mauguin symbol in short notation (e.g., "Fd-3m" for space group 227).
    pub hm_symbol: *const c_char,
    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Symmetry operations in the input cell.
    pub operations: MoyoOperations,
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// Spglib's `crystallographic_orbits` not `equivalent_atoms`
    /// The `i`th atom in the input cell is equivalent to the `orbits[i]`th atom in the **input** cell.
    /// For example, orbits=[0, 0, 2, 2, 2, 2] means the first two atoms are equivalent and the last four atoms are equivalent to each other.
    pub orbits: *const usize,
    /// Wyckoff letters for each site in the input cell.
    pub wyckoffs: *const c_char,
    /// Site symmetry symbols for each site in the input cell.
    /// The orientation of the site symmetry is w.r.t. the standardized cell.
    pub site_symmetry_symbols: *const *const c_char,
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
        // Identification
        let hm_symbol = CString::new(dataset.hm_symbol).expect("CString::new failed");

        // Site symmetry
        let wyckoffs = CString::new(dataset.wyckoffs.into_iter().collect::<String>())
            .expect("CString::new failed");
        let site_symmetry_symbols_cstring = dataset
            .site_symmetry_symbols
            .iter()
            .map(|s| CString::new(s.as_str()).expect("CString::new failed"))
            .collect::<Vec<_>>();
        let mut site_symmetry_symbols_ptr = Vec::with_capacity(site_symmetry_symbols_cstring.len());
        for s in site_symmetry_symbols_cstring {
            site_symmetry_symbols_ptr.push(s.into_raw() as *const c_char);
        }

        Self {
            // Identification
            number: dataset.number,
            hall_number: dataset.hall_number,
            hm_symbol: hm_symbol.into_raw(),
            // Symmetry operations in the input cell
            operations: (&dataset.operations).into(),
            // Site symmetry
            orbits: dataset.orbits.leak().as_ptr(),
            wyckoffs: wyckoffs.into_raw(),
            site_symmetry_symbols: site_symmetry_symbols_ptr.leak().as_ptr(),
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
pub extern "C" fn free_moyo_dataset(dataset: *mut MoyoDataset) {
    if dataset.is_null() {
        return;
    }
    unsafe {
        let dataset = Box::from_raw(dataset);
        if !dataset.hm_symbol.is_null() {
            drop(CString::from_raw(dataset.hm_symbol as *mut c_char));
        }
        free_moyo_operations(dataset.operations);
    }
}
