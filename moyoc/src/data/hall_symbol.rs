use core::ffi::c_char;

use moyo::data::{hall_symbol_entry, layer_hall_symbol_entry};

use crate::ffi::{free_cstring, leak_cstring};

/// An entry containing space-group information for a specified Hall number,
/// created by `moyo_hall_symbol_entry_new` and freed by
/// `moyo_hall_symbol_entry_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoHallSymbolEntry {
    /// Number for Hall symbols (1 - 530).
    pub hall_number: i32,
    /// ITA number for space group types (1 - 230).
    pub number: i32,
    /// Number for arithmetic crystal classes (1 - 73).
    pub arithmetic_number: i32,
    /// Setting.
    pub setting: *const c_char,
    /// Hall symbol.
    pub hall_symbol: *const c_char,
    /// Hermann-Mauguin symbol in short notation.
    pub hm_short: *const c_char,
    /// Hermann-Mauguin symbol in full notation.
    pub hm_full: *const c_char,
    /// Centering symbol ("P", "A", "B", "C", "I", "R", or "F").
    pub centering: *const c_char,
}

/// Look up the space-group Hall symbol entry for `hall_number` (1 - 530).
/// Returns NULL if `hall_number` is out of range.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_hall_symbol_entry_new(hall_number: i32) -> *mut MoyoHallSymbolEntry {
    // Reject lower-bound invalid numbers here: moyo indexes the database with
    // `hall_number - 1`, which overflows for extreme negative values.
    if hall_number < 1 {
        return std::ptr::null_mut();
    }
    match hall_symbol_entry(hall_number) {
        Some(entry) => Box::into_raw(Box::new(MoyoHallSymbolEntry {
            hall_number: entry.hall_number,
            number: entry.number,
            arithmetic_number: entry.arithmetic_number,
            setting: leak_cstring(entry.setting.to_string()),
            hall_symbol: leak_cstring(entry.hall_symbol.to_string()),
            hm_short: leak_cstring(entry.hm_short.to_string()),
            hm_full: leak_cstring(entry.hm_full.to_string()),
            centering: leak_cstring(format!("{:?}", entry.centering)),
        })),
        None => std::ptr::null_mut(),
    }
}

/// Free an entry created by `moyo_hall_symbol_entry_new`. Passing NULL is a no-op.
///
/// # Safety
/// `entry` must be a pointer returned by `moyo_hall_symbol_entry_new` that has
/// not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_hall_symbol_entry_free(entry: *mut MoyoHallSymbolEntry) {
    if entry.is_null() {
        return;
    }
    unsafe {
        let entry = Box::from_raw(entry);
        free_cstring(entry.setting);
        free_cstring(entry.hall_symbol);
        free_cstring(entry.hm_short);
        free_cstring(entry.hm_full);
        free_cstring(entry.centering);
    }
}

/// An entry containing layer-group information for a specified layer Hall
/// number, created by `moyo_layer_hall_symbol_entry_new` and freed by
/// `moyo_layer_hall_symbol_entry_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoLayerHallSymbolEntry {
    /// Number for layer Hall symbols (1 - 116).
    pub hall_number: i32,
    /// Layer-group number (1 - 80).
    pub number: i32,
    /// Number for layer arithmetic crystal classes (1 - 43).
    pub arithmetic_number: i32,
    /// Setting.
    pub setting: *const c_char,
    /// Hall symbol with lowercase `p`/`c` lattice prefix (layer convention).
    pub hall_symbol: *const c_char,
    /// Hermann-Mauguin symbol in short notation.
    pub hm_short: *const c_char,
    /// Hermann-Mauguin symbol in full notation.
    pub hm_full: *const c_char,
    /// Layer centering symbol ("P" or "C").
    pub centering: *const c_char,
}

/// Look up the layer-group Hall symbol entry for `hall_number` (1 - 116).
/// Returns NULL if `hall_number` is out of range.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_layer_hall_symbol_entry_new(
    hall_number: i32,
) -> *mut MoyoLayerHallSymbolEntry {
    if hall_number < 1 {
        return std::ptr::null_mut();
    }
    match layer_hall_symbol_entry(hall_number) {
        Some(entry) => Box::into_raw(Box::new(MoyoLayerHallSymbolEntry {
            hall_number: entry.hall_number,
            number: entry.number,
            arithmetic_number: entry.arithmetic_number,
            setting: leak_cstring(entry.setting.to_string()),
            hall_symbol: leak_cstring(entry.hall_symbol.to_string()),
            hm_short: leak_cstring(entry.hm_short.to_string()),
            hm_full: leak_cstring(entry.hm_full.to_string()),
            centering: leak_cstring(format!("{:?}", entry.centering)),
        })),
        None => std::ptr::null_mut(),
    }
}

/// Free an entry created by `moyo_layer_hall_symbol_entry_new`. Passing NULL is a no-op.
///
/// # Safety
/// `entry` must be a pointer returned by `moyo_layer_hall_symbol_entry_new`
/// that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_layer_hall_symbol_entry_free(entry: *mut MoyoLayerHallSymbolEntry) {
    if entry.is_null() {
        return;
    }
    unsafe {
        let entry = Box::from_raw(entry);
        free_cstring(entry.setting);
        free_cstring(entry.hall_symbol);
        free_cstring(entry.hm_short);
        free_cstring(entry.hm_full);
        free_cstring(entry.centering);
    }
}
