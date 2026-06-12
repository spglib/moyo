use moyo::base::{AngleTolerance, Lattice, LayerLattice, Operations};
use moyo::data::LayerSetting;
use moyo::identify::LayerGroup;

use crate::base::MoyoOperations;
use crate::data::MoyoLayerSetting;
use crate::identify::epsilon_or_default;

/// Layer group identified by `moyo_layer_group_new` and freed by
/// `moyo_layer_group_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoLayerGroup {
    /// Layer-group number for the identified layer group (1 - 80).
    pub number: i32,
    /// Layer Hall symbol number (1 - 116) for the chosen setting.
    pub hall_number: i32,
    /// Linear part of the transformation from the input primitive basis to the
    /// standardized basis.
    pub linear: [[i32; 3]; 3],
    /// Origin shift of the transformation from the input primitive basis to the
    /// standardized basis.
    pub origin_shift: [f64; 3],
}

/// Identify the layer group from primitive layer-cell rotations and translations.
///
/// - `prim_rotations`, `prim_translations`: rotation and translation parts of
///   the `num_operations` layer-group operations in the primitive layer cell.
/// - `basis`: row-wise basis vectors of the primitive layer lattice, `basis[i]`
///   is the `i`th basis vector. `a, b` must lie in the xy-plane and `c` along z
///   (the layer-group periodicity contract). Pass NULL to assume an identity
///   basis.
/// - `hall_number`: preference for layer Hall symbol, only used with
///   `MOYO_LAYER_SETTING_HALL_NUMBER`.
/// - `epsilon`: numerical tolerance for matching translations. Pass a negative
///   value to use the default tolerance (1e-4).
/// - Returns NULL if the identification fails or the arguments are invalid.
///   Free the result with `moyo_layer_group_free`.
///
/// # Safety
/// `prim_rotations` and `prim_translations` must point to `num_operations`
/// elements each, and `basis` must point to 3 basis vectors or be NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_layer_group_new(
    prim_rotations: *const [[i32; 3]; 3],
    prim_translations: *const [f64; 3],
    num_operations: i32,
    basis: *const [f64; 3],
    setting: MoyoLayerSetting,
    hall_number: i32,
    epsilon: f64,
) -> *mut MoyoLayerGroup {
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
        MoyoLayerSetting::HallNumber => {
            if hall_number <= 0 {
                return std::ptr::null_mut();
            }
            LayerSetting::HallNumber(hall_number)
        }
        MoyoLayerSetting::Spglib => LayerSetting::Spglib,
        MoyoLayerSetting::Standard => LayerSetting::Standard,
    };
    let epsilon = epsilon_or_default(epsilon);

    let layer_group = if basis.is_null() {
        LayerGroup::new(&prim_operations, setting, epsilon)
    } else {
        let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
        LayerLattice::new(Lattice::from_basis(basis), epsilon, AngleTolerance::Default).and_then(
            |layer_lattice| {
                LayerGroup::from_lattice(&layer_lattice, &prim_operations, setting, epsilon)
            },
        )
    };
    match layer_group {
        Ok(layer_group) => Box::into_raw(Box::new(MoyoLayerGroup {
            number: layer_group.number,
            hall_number: layer_group.hall_number,
            linear: layer_group.transformation.linear_as_array(),
            origin_shift: layer_group.transformation.origin_shift_as_array(),
        })),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a layer group created by `moyo_layer_group_new`. Passing NULL is a no-op.
///
/// # Safety
/// `layer_group` must be a pointer returned by `moyo_layer_group_new` that has
/// not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_layer_group_free(layer_group: *mut MoyoLayerGroup) {
    if layer_group.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(layer_group));
    }
}
