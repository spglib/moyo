use moyo::MoyoMagneticDataset as MagneticDataset;
use moyo::base::{
    AngleTolerance, Collinear, MagneticCell, NonCollinear, RotationMagneticMomentAction,
};
use moyo::utils::{to_3_slice, to_3x3_slice};

use crate::base::magnetic_cell::{
    moyo_collinear_magnetic_cell_free_members, moyo_noncollinear_magnetic_cell_free_members,
};
use crate::base::magnetic_operation::moyo_magnetic_operations_free_members;
use crate::base::{
    MoyoCell, MoyoCollinearMagneticCell, MoyoMagneticOperations, MoyoNonCollinearMagneticCell,
};
use crate::ffi::{free_slice, leak_slice};

/// A dataset of the magnetic symmetry analysis of a collinear magnetic
/// structure, created by `moyo_collinear_magnetic_dataset_new` and freed by
/// `moyo_collinear_magnetic_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoCollinearMagneticDataset {
    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------
    /// Serial number of UNI (and BNS) symbols.
    pub uni_number: i32,
    // ------------------------------------------------------------------------
    // Input magnetic cell
    // ------------------------------------------------------------------------
    /// Number of atoms in the input magnetic cell.
    /// Arrays `orbits` and `mapping_std_prim` have this length.
    pub num_atoms: i32,
    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Magnetic symmetry operations in the input cell.
    pub magnetic_operations: MoyoMagneticOperations,
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
    /// in the **input** magnetic cell.
    pub orbits: *const i32,
    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Standardized magnetic cell
    pub std_mag_cell: MoyoCollinearMagneticCell,
    /// Linear part of transformation from the input magnetic cell to the standardized one.
    pub std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input magnetic cell to the standardized one.
    pub std_origin_shift: [f64; 3],
    /// Rigid rotation
    pub std_rotation_matrix: [[f64; 3]; 3],
    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Primitive standardized magnetic cell
    pub prim_std_mag_cell: MoyoCollinearMagneticCell,
    /// Linear part of transformation from the input magnetic cell to the primitive standardized one.
    pub prim_std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input magnetic cell to the primitive standardized one.
    pub prim_std_origin_shift: [f64; 3],
    /// Mapping sites in the input magnetic cell to those in the primitive standardized one.
    /// The `i`th atom in the input magnetic cell is mapped to the `mapping_std_prim[i]`th
    /// atom in the primitive standardized magnetic cell.
    pub mapping_std_prim: *const i32,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    /// -1 if the default angle tolerance is used.
    pub angle_tolerance: f64,
    /// Actually used `mag_symprec` in iterative symmetry search.
    pub mag_symprec: f64,
}

impl MoyoCollinearMagneticDataset {
    fn new(dataset: MagneticDataset<Collinear>, num_atoms: i32) -> Self {
        let orbits = dataset.orbits.iter().map(|&x| x as i32).collect();
        let mapping_std_prim = dataset.mapping_std_prim.iter().map(|&x| x as i32).collect();

        Self {
            // Magnetic space-group type
            uni_number: dataset.uni_number,
            // Input magnetic cell
            num_atoms,
            // Magnetic symmetry operations in the input cell
            magnetic_operations: (&dataset.magnetic_operations).into(),
            // Site symmetry
            orbits: leak_slice(orbits),
            // Standardized magnetic cell
            std_mag_cell: (&dataset.std_mag_cell).into(),
            std_linear: to_3x3_slice(&dataset.std_linear),
            std_origin_shift: to_3_slice(&dataset.std_origin_shift),
            std_rotation_matrix: to_3x3_slice(&dataset.std_rotation_matrix),
            // Primitive standardized magnetic cell
            prim_std_mag_cell: (&dataset.prim_std_mag_cell).into(),
            prim_std_linear: to_3x3_slice(&dataset.prim_std_linear),
            prim_std_origin_shift: to_3_slice(&dataset.prim_std_origin_shift),
            mapping_std_prim: leak_slice(mapping_std_prim),
            // Final parameters
            symprec: dataset.symprec,
            angle_tolerance: angle_tolerance_to_f64(dataset.angle_tolerance),
            mag_symprec: dataset.mag_symprec,
        }
    }
}

/// A dataset of the magnetic symmetry analysis of a non-collinear magnetic
/// structure, created by `moyo_noncollinear_magnetic_dataset_new` and freed by
/// `moyo_noncollinear_magnetic_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoNonCollinearMagneticDataset {
    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------
    /// Serial number of UNI (and BNS) symbols.
    pub uni_number: i32,
    // ------------------------------------------------------------------------
    // Input magnetic cell
    // ------------------------------------------------------------------------
    /// Number of atoms in the input magnetic cell.
    /// Arrays `orbits` and `mapping_std_prim` have this length.
    pub num_atoms: i32,
    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Magnetic symmetry operations in the input cell.
    pub magnetic_operations: MoyoMagneticOperations,
    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// The `i`th atom in the input magnetic cell is equivalent to the `orbits[i]`th atom
    /// in the **input** magnetic cell.
    pub orbits: *const i32,
    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Standardized magnetic cell
    pub std_mag_cell: MoyoNonCollinearMagneticCell,
    /// Linear part of transformation from the input magnetic cell to the standardized one.
    pub std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input magnetic cell to the standardized one.
    pub std_origin_shift: [f64; 3],
    /// Rigid rotation
    pub std_rotation_matrix: [[f64; 3]; 3],
    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Primitive standardized magnetic cell
    pub prim_std_mag_cell: MoyoNonCollinearMagneticCell,
    /// Linear part of transformation from the input magnetic cell to the primitive standardized one.
    pub prim_std_linear: [[f64; 3]; 3],
    /// Origin shift of transformation from the input magnetic cell to the primitive standardized one.
    pub prim_std_origin_shift: [f64; 3],
    /// Mapping sites in the input magnetic cell to those in the primitive standardized one.
    /// The `i`th atom in the input magnetic cell is mapped to the `mapping_std_prim[i]`th
    /// atom in the primitive standardized magnetic cell.
    pub mapping_std_prim: *const i32,
    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actually used `symprec` in iterative symmetry search.
    pub symprec: f64,
    /// Actually used `angle_tolerance` in iterative symmetry search.
    /// -1 if the default angle tolerance is used.
    pub angle_tolerance: f64,
    /// Actually used `mag_symprec` in iterative symmetry search.
    pub mag_symprec: f64,
}

impl MoyoNonCollinearMagneticDataset {
    fn new(dataset: MagneticDataset<NonCollinear>, num_atoms: i32) -> Self {
        let orbits = dataset.orbits.iter().map(|&x| x as i32).collect();
        let mapping_std_prim = dataset.mapping_std_prim.iter().map(|&x| x as i32).collect();

        Self {
            // Magnetic space-group type
            uni_number: dataset.uni_number,
            // Input magnetic cell
            num_atoms,
            // Magnetic symmetry operations in the input cell
            magnetic_operations: (&dataset.magnetic_operations).into(),
            // Site symmetry
            orbits: leak_slice(orbits),
            // Standardized magnetic cell
            std_mag_cell: (&dataset.std_mag_cell).into(),
            std_linear: to_3x3_slice(&dataset.std_linear),
            std_origin_shift: to_3_slice(&dataset.std_origin_shift),
            std_rotation_matrix: to_3x3_slice(&dataset.std_rotation_matrix),
            // Primitive standardized magnetic cell
            prim_std_mag_cell: (&dataset.prim_std_mag_cell).into(),
            prim_std_linear: to_3x3_slice(&dataset.prim_std_linear),
            prim_std_origin_shift: to_3_slice(&dataset.prim_std_origin_shift),
            mapping_std_prim: leak_slice(mapping_std_prim),
            // Final parameters
            symprec: dataset.symprec,
            angle_tolerance: angle_tolerance_to_f64(dataset.angle_tolerance),
            mag_symprec: dataset.mag_symprec,
        }
    }
}

/// Analyze the magnetic symmetry of the given collinear magnetic cell and create a dataset.
///
/// - `basis`: row-wise basis vectors, `basis[i]` is the `i`th basis vector.
/// - `positions`: fractional coordinates of the `num_atoms` sites.
/// - `numbers`: atomic numbers of the `num_atoms` sites.
/// - `magnetic_moments`: collinear magnetic moments of the `num_atoms` sites.
/// - `angle_tolerance`: tolerance of angle between basis vectors in radians.
///   Pass a negative value to use the default angle tolerance.
/// - `mag_symprec`: tolerance for matching magnetic moments.
///   Pass a negative value to reuse `symprec`.
/// - `is_axial`: whether magnetic moments transform as axial vectors.
/// - Returns NULL if the symmetry search fails or the arguments are invalid.
///
/// # Safety
/// `basis` must point to 3 basis vectors, and `positions`, `numbers`, and
/// `magnetic_moments` must point to `num_atoms` elements each.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_collinear_magnetic_dataset_new(
    basis: *const [f64; 3],
    positions: *const [f64; 3],
    numbers: *const i32,
    magnetic_moments: *const f64,
    num_atoms: i32,
    symprec: f64,
    angle_tolerance: f64,
    mag_symprec: f64,
    is_axial: bool,
    rotate_basis: bool,
) -> *mut MoyoCollinearMagneticDataset {
    if basis.is_null()
        || positions.is_null()
        || numbers.is_null()
        || magnetic_moments.is_null()
        || num_atoms <= 0
    {
        return std::ptr::null_mut();
    }
    let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
    let magnetic_cell = MoyoCollinearMagneticCell {
        cell: MoyoCell {
            basis,
            positions,
            numbers,
            num_atoms,
        },
        magnetic_moments,
    };
    let magnetic_cell: MagneticCell<Collinear> = (&magnetic_cell).into();

    let dataset = MagneticDataset::new(
        &magnetic_cell,
        symprec,
        angle_tolerance_from_f64(angle_tolerance),
        mag_symprec_from_f64(mag_symprec),
        action_from_is_axial(is_axial),
        rotate_basis,
    );
    match dataset {
        Ok(dataset) => Box::into_raw(Box::new(MoyoCollinearMagneticDataset::new(
            dataset, num_atoms,
        ))),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a dataset created by `moyo_collinear_magnetic_dataset_new`. Passing NULL is a no-op.
///
/// # Safety
/// `dataset` must be a pointer returned by `moyo_collinear_magnetic_dataset_new`
/// that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_collinear_magnetic_dataset_free(
    dataset: *mut MoyoCollinearMagneticDataset,
) {
    if dataset.is_null() {
        return;
    }
    unsafe {
        let dataset = Box::from_raw(dataset);
        let num_atoms = dataset.num_atoms as usize;

        // Magnetic symmetry operations in the input cell
        moyo_magnetic_operations_free_members(&dataset.magnetic_operations);

        // Site symmetry
        free_slice(dataset.orbits, num_atoms);

        // Standardized magnetic cell
        moyo_collinear_magnetic_cell_free_members(&dataset.std_mag_cell);

        // Primitive standardized magnetic cell
        moyo_collinear_magnetic_cell_free_members(&dataset.prim_std_mag_cell);
        free_slice(dataset.mapping_std_prim, num_atoms);
    }
}

/// Analyze the magnetic symmetry of the given non-collinear magnetic cell and create a dataset.
///
/// - `basis`: row-wise basis vectors, `basis[i]` is the `i`th basis vector.
/// - `positions`: fractional coordinates of the `num_atoms` sites.
/// - `numbers`: atomic numbers of the `num_atoms` sites.
/// - `magnetic_moments`: non-collinear magnetic moments of the `num_atoms` sites
///   in Cartesian coordinates.
/// - `angle_tolerance`: tolerance of angle between basis vectors in radians.
///   Pass a negative value to use the default angle tolerance.
/// - `mag_symprec`: tolerance for matching magnetic moments.
///   Pass a negative value to reuse `symprec`.
/// - `is_axial`: whether magnetic moments transform as axial vectors.
/// - Returns NULL if the symmetry search fails or the arguments are invalid.
///
/// # Safety
/// `basis` must point to 3 basis vectors, and `positions`, `numbers`, and
/// `magnetic_moments` must point to `num_atoms` elements each.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_noncollinear_magnetic_dataset_new(
    basis: *const [f64; 3],
    positions: *const [f64; 3],
    numbers: *const i32,
    magnetic_moments: *const [f64; 3],
    num_atoms: i32,
    symprec: f64,
    angle_tolerance: f64,
    mag_symprec: f64,
    is_axial: bool,
    rotate_basis: bool,
) -> *mut MoyoNonCollinearMagneticDataset {
    if basis.is_null()
        || positions.is_null()
        || numbers.is_null()
        || magnetic_moments.is_null()
        || num_atoms <= 0
    {
        return std::ptr::null_mut();
    }
    let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
    let magnetic_cell = MoyoNonCollinearMagneticCell {
        cell: MoyoCell {
            basis,
            positions,
            numbers,
            num_atoms,
        },
        magnetic_moments,
    };
    let magnetic_cell: MagneticCell<NonCollinear> = (&magnetic_cell).into();

    let dataset = MagneticDataset::new(
        &magnetic_cell,
        symprec,
        angle_tolerance_from_f64(angle_tolerance),
        mag_symprec_from_f64(mag_symprec),
        action_from_is_axial(is_axial),
        rotate_basis,
    );
    match dataset {
        Ok(dataset) => Box::into_raw(Box::new(MoyoNonCollinearMagneticDataset::new(
            dataset, num_atoms,
        ))),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a dataset created by `moyo_noncollinear_magnetic_dataset_new`. Passing NULL is a no-op.
///
/// # Safety
/// `dataset` must be a pointer returned by `moyo_noncollinear_magnetic_dataset_new`
/// that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_noncollinear_magnetic_dataset_free(
    dataset: *mut MoyoNonCollinearMagneticDataset,
) {
    if dataset.is_null() {
        return;
    }
    unsafe {
        let dataset = Box::from_raw(dataset);
        let num_atoms = dataset.num_atoms as usize;

        // Magnetic symmetry operations in the input cell
        moyo_magnetic_operations_free_members(&dataset.magnetic_operations);

        // Site symmetry
        free_slice(dataset.orbits, num_atoms);

        // Standardized magnetic cell
        moyo_noncollinear_magnetic_cell_free_members(&dataset.std_mag_cell);

        // Primitive standardized magnetic cell
        moyo_noncollinear_magnetic_cell_free_members(&dataset.prim_std_mag_cell);
        free_slice(dataset.mapping_std_prim, num_atoms);
    }
}

fn angle_tolerance_from_f64(angle_tolerance: f64) -> AngleTolerance {
    if angle_tolerance < 0.0 {
        AngleTolerance::default()
    } else {
        AngleTolerance::Radian(angle_tolerance)
    }
}

fn angle_tolerance_to_f64(angle_tolerance: AngleTolerance) -> f64 {
    if let AngleTolerance::Radian(angle_tolerance) = angle_tolerance {
        angle_tolerance
    } else {
        -1.0
    }
}

fn mag_symprec_from_f64(mag_symprec: f64) -> Option<f64> {
    if mag_symprec < 0.0 {
        None
    } else {
        Some(mag_symprec)
    }
}

fn action_from_is_axial(is_axial: bool) -> RotationMagneticMomentAction {
    if is_axial {
        RotationMagneticMomentAction::Axial
    } else {
        RotationMagneticMomentAction::Polar
    }
}
