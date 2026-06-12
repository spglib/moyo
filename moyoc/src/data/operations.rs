use moyo::data::{
    LayerSetting, Setting, magnetic_operations_from_uni_number, operations_from_layer_number,
    operations_from_number,
};

use crate::base::magnetic_operation::moyo_magnetic_operations_free_members;
use crate::base::operation::moyo_operations_free_members;
use crate::base::{MoyoMagneticOperations, MoyoOperations};
use crate::data::{MoyoLayerSetting, MoyoSetting};

/// Return symmetry operations for the given space-group ITA `number` (1 - 230).
///
/// - `hall_number`: preference for Hall symbol, only used with `MOYO_SETTING_HALL_NUMBER`
///   (`number` is ignored in that case).
/// - `primitive`: if true, return operations of the primitive cell; otherwise return
///   operations of the conventional cell.
/// - Returns NULL if the arguments are invalid. Free the returned operations with
///   `moyo_operations_free`.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_operations_from_number(
    number: i32,
    setting: MoyoSetting,
    hall_number: i32,
    primitive: bool,
) -> *mut MoyoOperations {
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
    match operations_from_number(number, setting, primitive) {
        Ok(operations) => Box::into_raw(Box::new((&operations).into())),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Return symmetry operations for the given layer-group `number` (1 - 80).
///
/// - `hall_number`: preference for layer Hall symbol, only used with
///   `MOYO_LAYER_SETTING_HALL_NUMBER` (`number` is ignored in that case).
/// - `primitive`: if true, return operations of the primitive cell; otherwise return
///   operations of the conventional cell.
/// - Returns NULL if the arguments are invalid. Free the returned operations with
///   `moyo_operations_free`.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_operations_from_layer_number(
    number: i32,
    setting: MoyoLayerSetting,
    hall_number: i32,
    primitive: bool,
) -> *mut MoyoOperations {
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
    match operations_from_layer_number(number, setting, primitive) {
        Ok(operations) => Box::into_raw(Box::new((&operations).into())),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free operations returned by `moyo_operations_from_number` or
/// `moyo_operations_from_layer_number`. Passing NULL is a no-op.
/// Do not call this on operations embedded in a dataset; those are freed by
/// the dataset's free function.
///
/// # Safety
/// `operations` must be a pointer returned by `moyo_operations_from_number` or
/// `moyo_operations_from_layer_number` that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_operations_free(operations: *mut MoyoOperations) {
    if operations.is_null() {
        return;
    }
    unsafe {
        let operations = Box::from_raw(operations);
        moyo_operations_free_members(&operations);
    }
}

/// Return magnetic symmetry operations for the given UNI (and BNS) serial
/// number `uni_number` (1 - 1651).
///
/// - `primitive`: if true, return magnetic operations of the primitive cell; otherwise
///   return operations of the conventional cell.
/// - Returns NULL if the arguments are invalid. Free the returned operations with
///   `moyo_magnetic_operations_free`.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_magnetic_operations_from_uni_number(
    uni_number: i32,
    primitive: bool,
) -> *mut MoyoMagneticOperations {
    match magnetic_operations_from_uni_number(uni_number, primitive) {
        Ok(magnetic_operations) => Box::into_raw(Box::new((&magnetic_operations).into())),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free magnetic operations returned by `moyo_magnetic_operations_from_uni_number`.
/// Passing NULL is a no-op.
/// Do not call this on magnetic operations embedded in a dataset; those are
/// freed by the dataset's free function.
///
/// # Safety
/// `magnetic_operations` must be a pointer returned by
/// `moyo_magnetic_operations_from_uni_number` that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_magnetic_operations_free(
    magnetic_operations: *mut MoyoMagneticOperations,
) {
    if magnetic_operations.is_null() {
        return;
    }
    unsafe {
        let magnetic_operations = Box::from_raw(magnetic_operations);
        moyo_magnetic_operations_free_members(&magnetic_operations);
    }
}
