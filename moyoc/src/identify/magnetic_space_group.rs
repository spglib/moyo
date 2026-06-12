use moyo::base::{Lattice, MagneticOperations};
use moyo::identify::MagneticSpaceGroup;

use crate::base::MoyoMagneticOperations;
use crate::identify::epsilon_or_default;

/// Magnetic space group identified by `moyo_magnetic_space_group_new` and
/// freed by `moyo_magnetic_space_group_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoMagneticSpaceGroup {
    /// Serial number of UNI (and BNS) symbols (1 - 1651).
    pub uni_number: i32,
    /// Linear part of the transformation from the input primitive basis to the
    /// standardized basis.
    pub linear: [[i32; 3]; 3],
    /// Origin shift of the transformation from the input primitive basis to the
    /// standardized basis.
    pub origin_shift: [f64; 3],
}

/// Identify the magnetic space group from primitive magnetic operations.
///
/// - `prim_rotations`, `prim_translations`, `prim_time_reversals`: rotation,
///   translation, and time-reversal parts of the `num_operations` magnetic
///   operations in the primitive cell.
/// - `basis`: row-wise basis vectors of the primitive lattice, `basis[i]` is
///   the `i`th basis vector. Pass NULL to assume an identity basis.
/// - `epsilon`: numerical tolerance for matching translations. Pass a negative
///   value to use the default tolerance (1e-4).
/// - Returns NULL if the identification fails or the arguments are invalid.
///   Free the result with `moyo_magnetic_space_group_free`.
///
/// # Safety
/// `prim_rotations`, `prim_translations`, and `prim_time_reversals` must point
/// to `num_operations` elements each, and `basis` must point to 3 basis
/// vectors or be NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_magnetic_space_group_new(
    prim_rotations: *const [[i32; 3]; 3],
    prim_translations: *const [f64; 3],
    prim_time_reversals: *const bool,
    num_operations: i32,
    basis: *const [f64; 3],
    epsilon: f64,
) -> *mut MoyoMagneticSpaceGroup {
    if prim_rotations.is_null()
        || prim_translations.is_null()
        || prim_time_reversals.is_null()
        || num_operations <= 0
    {
        return std::ptr::null_mut();
    }
    let prim_mag_operations: MagneticOperations = (&MoyoMagneticOperations {
        rotations: prim_rotations,
        translations: prim_translations,
        time_reversals: prim_time_reversals,
        num_operations,
    })
        .into();
    let epsilon = epsilon_or_default(epsilon);

    let magnetic_space_group = if basis.is_null() {
        MagneticSpaceGroup::new(&prim_mag_operations, epsilon)
    } else {
        let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
        MagneticSpaceGroup::from_lattice(&Lattice::from_basis(basis), &prim_mag_operations, epsilon)
    };
    match magnetic_space_group {
        Ok(magnetic_space_group) => Box::into_raw(Box::new(MoyoMagneticSpaceGroup {
            uni_number: magnetic_space_group.uni_number,
            linear: magnetic_space_group.transformation.linear_as_array(),
            origin_shift: magnetic_space_group.transformation.origin_shift_as_array(),
        })),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a magnetic space group created by `moyo_magnetic_space_group_new`.
/// Passing NULL is a no-op.
///
/// # Safety
/// `magnetic_space_group` must be a pointer returned by
/// `moyo_magnetic_space_group_new` that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_magnetic_space_group_free(
    magnetic_space_group: *mut MoyoMagneticSpaceGroup,
) {
    if magnetic_space_group.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(magnetic_space_group));
    }
}
