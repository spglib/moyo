use core::ffi::c_char;

use moyo::data::{ConstructType, get_magnetic_space_group_type};

use crate::ffi::{free_cstring, leak_cstring};

/// Magnetic space-group type information, created by
/// `moyo_magnetic_space_group_type_new` and freed by
/// `moyo_magnetic_space_group_type_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoMagneticSpaceGroupType {
    /// Serial number of UNI (and BNS) symbols (1 - 1651).
    pub uni_number: i32,
    /// Serial number in Litvin's Magnetic group tables.
    pub litvin_number: i32,
    /// BNS number, e.g. "151.32".
    pub bns_number: *const c_char,
    /// OG number, e.g. "153.4.1270".
    pub og_number: *const c_char,
    /// ITA number for the reference space group in BNS setting.
    pub number: i32,
    /// Construct type of the magnetic space group, from 1 to 4.
    pub construct_type: i32,
}

/// Look up the magnetic space-group type information for `uni_number` (1 - 1651).
/// Returns NULL if `uni_number` is out of range.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_magnetic_space_group_type_new(
    uni_number: i32,
) -> *mut MoyoMagneticSpaceGroupType {
    match get_magnetic_space_group_type(uni_number) {
        Some(magnetic_space_group_type) => {
            let construct_type = match magnetic_space_group_type.construct_type {
                ConstructType::Type1 => 1,
                ConstructType::Type2 => 2,
                ConstructType::Type3 => 3,
                ConstructType::Type4 => 4,
            };
            Box::into_raw(Box::new(MoyoMagneticSpaceGroupType {
                uni_number: magnetic_space_group_type.uni_number,
                litvin_number: magnetic_space_group_type.litvin_number,
                bns_number: leak_cstring(magnetic_space_group_type.bns_number.to_string()),
                og_number: leak_cstring(magnetic_space_group_type.og_number.to_string()),
                number: magnetic_space_group_type.number,
                construct_type,
            }))
        }
        None => std::ptr::null_mut(),
    }
}

/// Free a magnetic space-group type created by `moyo_magnetic_space_group_type_new`.
/// Passing NULL is a no-op.
///
/// # Safety
/// `magnetic_space_group_type` must be a pointer returned by
/// `moyo_magnetic_space_group_type_new` that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_magnetic_space_group_type_free(
    magnetic_space_group_type: *mut MoyoMagneticSpaceGroupType,
) {
    if magnetic_space_group_type.is_null() {
        return;
    }
    unsafe {
        let magnetic_space_group_type = Box::from_raw(magnetic_space_group_type);
        free_cstring(magnetic_space_group_type.bns_number);
        free_cstring(magnetic_space_group_type.og_number);
    }
}
