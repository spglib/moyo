use core::ffi::c_char;
use std::ffi::CString;

use moyo::MoyoDataset as Dataset;
use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;
use moyo::utils::{to_3_slice, to_3x3_slice};

use crate::base::{MoyoCell, MoyoOperations, free_moyo_cell, free_moyo_operations};
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
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    /// Standardized cell
    pub std_cell: MoyoCell,
    /// Linear part of transformation from the input cell to the standardized cell.
    pub std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input cell to the standardized cell.
    pub std_origin_shift: [f64; 3],
    /// Rigid rotation
    pub std_rotation_matrix: [[f64; 3]; 3],
    /// Pearson symbol for standardized cell
    pub pearson_symbol: *const c_char,
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    /// Primitive standardized cell
    pub prim_std_cell: MoyoCell,
    /// Linear part of transformation from the input cell to the primitive standardized cell.
    pub prim_std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input cell to the primitive standardized cell.
    pub prim_std_origin_shift: [f64; 3],
    /// Mapping sites in the input cell to those in the primitive standardized cell.
    /// The `i`th atom in the input cell is mapped to the `mapping_to_std_prim[i]`th atom in the primitive standardized cell.
    pub mapping_std_prim: *const usize,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    pub angle_tolerance: f64,
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

        // Standardized cell
        let pearson_symbol = CString::new(dataset.pearson_symbol).expect("CString::new failed");

        // Final parameters
        let angle_tolerance =
            if let AngleTolerance::Radian(angle_tolerance) = dataset.angle_tolerance {
                angle_tolerance
            } else {
                -1.0
            };

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
            // Standardized cell
            std_cell: (&dataset.std_cell).into(),
            std_linear: to_3x3_slice(&dataset.std_linear),
            std_origin_shift: to_3_slice(&dataset.std_origin_shift),
            std_rotation_matrix: to_3x3_slice(&dataset.std_rotation_matrix),
            pearson_symbol: pearson_symbol.into_raw(),
            // Primitive standardized cell
            prim_std_cell: (&dataset.prim_std_cell).into(),
            prim_std_linear: to_3x3_slice(&dataset.prim_std_linear),
            prim_std_origin_shift: to_3_slice(&dataset.prim_std_origin_shift),
            mapping_std_prim: dataset.mapping_std_prim.leak().as_ptr(),
            // Final parameters
            symprec: dataset.symprec,
            angle_tolerance,
        }
    }
}

#[unsafe(no_mangle)]
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
    // TODO: refactor conversion between MoyoSetting and Setting
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
    match dataset {
        Ok(dataset) => {
            let moyoc_dataset: MoyoDataset = dataset.into();
            Box::into_raw(Box::new(moyoc_dataset))
        }
        Err(_) => std::ptr::null_mut(),
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn free_moyo_dataset(dataset: *mut MoyoDataset) {
    if dataset.is_null() {
        return;
    }
    unsafe {
        let dataset = Box::from_raw(dataset);

        // Identification
        if !dataset.hm_symbol.is_null() {
            drop(CString::from_raw(dataset.hm_symbol as *mut c_char));
        }

        // Symmetry operations in the input cell
        free_moyo_operations(dataset.operations);

        // Site symmetry
        if !dataset.orbits.is_null() {
            let _ = Vec::from_raw_parts(
                dataset.orbits as *mut usize,
                dataset.std_cell.num_atoms as usize,
                dataset.std_cell.num_atoms as usize,
            );
        }
        if !dataset.wyckoffs.is_null() {
            drop(CString::from_raw(dataset.wyckoffs as *mut c_char));
        }
        if !dataset.site_symmetry_symbols.is_null() {
            let site_symmetry_symbols_ptr = std::slice::from_raw_parts(
                dataset.site_symmetry_symbols,
                dataset.std_cell.num_atoms as usize,
            );
            for &s in site_symmetry_symbols_ptr {
                if !s.is_null() {
                    drop(CString::from_raw(s as *mut c_char));
                }
            }
            let _ = Vec::from_raw_parts(
                dataset.site_symmetry_symbols as *mut *const c_char,
                dataset.std_cell.num_atoms as usize,
                dataset.std_cell.num_atoms as usize,
            );
        }

        // Standardized cell
        free_moyo_cell(dataset.std_cell);
        if !dataset.pearson_symbol.is_null() {
            drop(CString::from_raw(dataset.pearson_symbol as *mut c_char));
        }

        // Primitive standardized cell
        if !dataset.mapping_std_prim.is_null() {
            let _ = Vec::from_raw_parts(
                dataset.mapping_std_prim as *mut usize,
                dataset.prim_std_cell.num_atoms as usize,
                dataset.prim_std_cell.num_atoms as usize,
            );
        }
        free_moyo_cell(dataset.prim_std_cell);
    }
}
