use core::ffi::c_char;

use moyo::MoyoDataset as Dataset;
use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;
use moyo::utils::{to_3_slice, to_3x3_slice};

use crate::base::cell::moyo_cell_free_members;
use crate::base::operation::moyo_operations_free_members;
use crate::base::{MoyoCell, MoyoOperations};
use crate::data::MoyoSetting;
use crate::ffi::{free_cstring, free_slice, leak_cstring, leak_slice};

/// A dataset of the symmetry analysis, created by `moyo_dataset_new` and freed
/// by `moyo_dataset_free`.
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
    // Input cell
    // ------------------------------------------------------------------------
    /// Number of atoms in the input cell.
    /// Arrays `orbits`, `wyckoffs`, `site_symmetry_symbols`, and `mapping_std_prim`
    /// have this length.
    pub num_atoms: i32,
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
    pub orbits: *const i32,
    /// Wyckoff letters for each site in the input cell as a NUL-terminated string.
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
    /// The `i`th atom in the input cell is mapped to the `mapping_std_prim[i]`th atom in the primitive standardized cell.
    pub mapping_std_prim: *const i32,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    /// -1 if the default angle tolerance is used.
    pub angle_tolerance: f64,
}

impl MoyoDataset {
    fn new(dataset: Dataset, num_atoms: i32) -> Self {
        // Site symmetry
        let wyckoffs = dataset.wyckoffs.into_iter().collect::<String>();
        let site_symmetry_symbols = dataset
            .site_symmetry_symbols
            .iter()
            .map(|s| leak_cstring(s.clone()))
            .collect::<Vec<_>>();
        let orbits = dataset.orbits.iter().map(|&x| x as i32).collect();
        let mapping_std_prim = dataset.mapping_std_prim.iter().map(|&x| x as i32).collect();

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
            hm_symbol: leak_cstring(dataset.hm_symbol),
            // Input cell
            num_atoms,
            // Symmetry operations in the input cell
            operations: (&dataset.operations).into(),
            // Site symmetry
            orbits: leak_slice(orbits),
            wyckoffs: leak_cstring(wyckoffs),
            site_symmetry_symbols: leak_slice(site_symmetry_symbols),
            // Standardized cell
            std_cell: (&dataset.std_cell).into(),
            std_linear: to_3x3_slice(&dataset.std_linear),
            std_origin_shift: to_3_slice(&dataset.std_origin_shift),
            std_rotation_matrix: to_3x3_slice(&dataset.std_rotation_matrix),
            pearson_symbol: leak_cstring(dataset.pearson_symbol),
            // Primitive standardized cell
            prim_std_cell: (&dataset.prim_std_cell).into(),
            prim_std_linear: to_3x3_slice(&dataset.prim_std_linear),
            prim_std_origin_shift: to_3_slice(&dataset.prim_std_origin_shift),
            mapping_std_prim: leak_slice(mapping_std_prim),
            // Final parameters
            symprec: dataset.symprec,
            angle_tolerance,
        }
    }
}

/// Analyze the symmetry of the given cell and create a dataset.
///
/// - `basis`: row-wise basis vectors, `basis[i]` is the `i`th basis vector.
/// - `positions`: fractional coordinates of the `num_atoms` sites.
/// - `numbers`: atomic numbers of the `num_atoms` sites.
/// - `angle_tolerance`: tolerance of angle between basis vectors in radians.
///   Pass a negative value to use the default angle tolerance.
/// - `hall_number`: preference for Hall symbol, only used with `MOYO_SETTING_HALL_NUMBER`.
/// - Returns NULL if the symmetry search fails or the arguments are invalid.
///
/// # Safety
/// `basis` must point to 3 basis vectors, and `positions` and `numbers` must
/// point to `num_atoms` elements each.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_dataset_new(
    basis: *const [f64; 3],
    positions: *const [f64; 3],
    numbers: *const i32,
    num_atoms: i32,
    symprec: f64,
    angle_tolerance: f64,
    setting: MoyoSetting,
    hall_number: i32,
    rotate_basis: bool,
) -> *mut MoyoDataset {
    if basis.is_null() || positions.is_null() || numbers.is_null() || num_atoms <= 0 {
        return std::ptr::null_mut();
    }
    let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
    let cell = MoyoCell {
        basis,
        positions,
        numbers,
        num_atoms,
    };
    let cell: Cell = (&cell).into();

    let angle_tolerance = if angle_tolerance < 0.0 {
        AngleTolerance::default()
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

    let dataset = Dataset::new(&cell, symprec, angle_tolerance, setting, rotate_basis);
    match dataset {
        Ok(dataset) => Box::into_raw(Box::new(MoyoDataset::new(dataset, num_atoms))),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a dataset created by `moyo_dataset_new`. Passing NULL is a no-op.
///
/// # Safety
/// `dataset` must be a pointer returned by `moyo_dataset_new` that has not
/// been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_dataset_free(dataset: *mut MoyoDataset) {
    if dataset.is_null() {
        return;
    }
    unsafe {
        let dataset = Box::from_raw(dataset);
        let num_atoms = dataset.num_atoms as usize;

        // Identification
        free_cstring(dataset.hm_symbol);

        // Symmetry operations in the input cell
        moyo_operations_free_members(&dataset.operations);

        // Site symmetry
        free_slice(dataset.orbits, num_atoms);
        free_cstring(dataset.wyckoffs);
        if !dataset.site_symmetry_symbols.is_null() {
            let site_symmetry_symbols =
                std::slice::from_raw_parts(dataset.site_symmetry_symbols, num_atoms);
            for &s in site_symmetry_symbols {
                free_cstring(s);
            }
            free_slice(dataset.site_symmetry_symbols, num_atoms);
        }

        // Standardized cell
        moyo_cell_free_members(&dataset.std_cell);
        free_cstring(dataset.pearson_symbol);

        // Primitive standardized cell
        moyo_cell_free_members(&dataset.prim_std_cell);
        free_slice(dataset.mapping_std_prim, num_atoms);
    }
}
