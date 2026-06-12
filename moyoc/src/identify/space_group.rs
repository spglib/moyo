use moyo::base::{Lattice, Operations};
use moyo::data::Setting;
use moyo::identify::SpaceGroup;

use crate::base::MoyoOperations;
use crate::data::MoyoSetting;
use crate::identify::epsilon_or_default;

/// Space group identified by `moyo_space_group_new` and freed by
/// `moyo_space_group_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoSpaceGroup {
    /// ITA number for the identified space group (1 - 230).
    pub number: i32,
    /// Hall symbol number (1 - 530) for the chosen setting.
    pub hall_number: i32,
    /// Linear part of the transformation from the input primitive basis to the
    /// standardized basis.
    pub linear: [[i32; 3]; 3],
    /// Origin shift of the transformation from the input primitive basis to the
    /// standardized basis.
    pub origin_shift: [f64; 3],
}

/// Identify the space group from primitive rotations and translations.
///
/// - `prim_rotations`, `prim_translations`: rotation and translation parts of
///   the `num_operations` operations in the primitive cell.
/// - `basis`: row-wise basis vectors of the primitive lattice, `basis[i]` is
///   the `i`th basis vector. Pass NULL to assume an identity basis.
/// - `hall_number`: preference for Hall symbol, only used with
///   `MOYO_SETTING_HALL_NUMBER`.
/// - `epsilon`: numerical tolerance for matching translations. Pass a negative
///   value to use the default tolerance (1e-4).
/// - Returns NULL if the identification fails or the arguments are invalid.
///   Free the result with `moyo_space_group_free`.
///
/// # Safety
/// `prim_rotations` and `prim_translations` must point to `num_operations`
/// elements each, and `basis` must point to 3 basis vectors or be NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_space_group_new(
    prim_rotations: *const [[i32; 3]; 3],
    prim_translations: *const [f64; 3],
    num_operations: i32,
    basis: *const [f64; 3],
    setting: MoyoSetting,
    hall_number: i32,
    epsilon: f64,
) -> *mut MoyoSpaceGroup {
    if prim_rotations.is_null() || prim_translations.is_null() || num_operations <= 0 {
        return std::ptr::null_mut();
    }
    let prim_operations: Operations = (&MoyoOperations {
        rotations: prim_rotations,
        translations: prim_translations,
        num_operations,
    })
        .into();
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
    let epsilon = epsilon_or_default(epsilon);

    let space_group = if basis.is_null() {
        SpaceGroup::new(&prim_operations, setting, epsilon)
    } else {
        let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
        SpaceGroup::from_lattice(
            &Lattice::from_basis(basis),
            &prim_operations,
            setting,
            epsilon,
        )
    };
    match space_group {
        Ok(space_group) => Box::into_raw(Box::new(MoyoSpaceGroup {
            number: space_group.number,
            hall_number: space_group.hall_number,
            linear: space_group.transformation.linear_as_array(),
            origin_shift: space_group.transformation.origin_shift_as_array(),
        })),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a space group created by `moyo_space_group_new`. Passing NULL is a no-op.
///
/// # Safety
/// `space_group` must be a pointer returned by `moyo_space_group_new` that has
/// not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_space_group_free(space_group: *mut MoyoSpaceGroup) {
    if space_group.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(space_group));
    }
}
