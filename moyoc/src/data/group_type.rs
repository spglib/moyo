use core::ffi::c_char;

use moyo::data::{
    CrystalFamily, CrystalSystem, LatticeSystem, LayerSetting, Setting,
    arithmetic_crystal_class_entry, hall_symbol_entry, layer_arithmetic_crystal_class_entry,
    layer_hall_symbol_entry,
};

use crate::ffi::{free_cstring, leak_cstring};

/// Space-group type information, created by `moyo_space_group_type_new` and
/// freed by `moyo_space_group_type_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoSpaceGroupType {
    /// ITA number for space group types (1 - 230).
    pub number: i32,
    /// Hermann-Mauguin symbol in short notation.
    pub hm_short: *const c_char,
    /// Hermann-Mauguin symbol in full notation.
    pub hm_full: *const c_char,
    /// Number for arithmetic crystal classes (1 - 73).
    pub arithmetic_number: i32,
    /// Symbol for arithmetic crystal class.
    pub arithmetic_symbol: *const c_char,
    /// Geometric crystal class.
    pub geometric_crystal_class: *const c_char,
    /// Crystal system.
    pub crystal_system: *const c_char,
    /// Bravais class.
    pub bravais_class: *const c_char,
    /// Lattice system.
    pub lattice_system: *const c_char,
    /// Crystal family.
    pub crystal_family: *const c_char,
}

/// Look up the space-group type information for ITA `number` (1 - 230).
/// Returns NULL if `number` is out of range.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_space_group_type_new(number: i32) -> *mut MoyoSpaceGroupType {
    let Some(hall_number) = Setting::default().hall_number(number) else {
        return std::ptr::null_mut();
    };
    let hall_symbol = hall_symbol_entry(hall_number).unwrap();
    let arithmetic_entry = arithmetic_crystal_class_entry(hall_symbol.arithmetic_number).unwrap();

    let geometric_crystal_class = arithmetic_entry.geometric_crystal_class;
    let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);
    let bravais_class = arithmetic_entry.bravais_class;
    let lattice_system = LatticeSystem::from_bravais_class(bravais_class);
    let crystal_family = CrystalFamily::from_lattice_system(lattice_system);

    Box::into_raw(Box::new(MoyoSpaceGroupType {
        number,
        hm_short: leak_cstring(hall_symbol.hm_short.to_string()),
        hm_full: leak_cstring(hall_symbol.hm_full.to_string()),
        arithmetic_number: hall_symbol.arithmetic_number,
        arithmetic_symbol: leak_cstring(arithmetic_entry.symbol.to_string()),
        geometric_crystal_class: leak_cstring(geometric_crystal_class.to_string()),
        crystal_system: leak_cstring(crystal_system.to_string()),
        bravais_class: leak_cstring(bravais_class.to_string()),
        lattice_system: leak_cstring(lattice_system.to_string()),
        crystal_family: leak_cstring(crystal_family.to_string()),
    }))
}

/// Free a space-group type created by `moyo_space_group_type_new`.
/// Passing NULL is a no-op.
///
/// # Safety
/// `space_group_type` must be a pointer returned by `moyo_space_group_type_new`
/// that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_space_group_type_free(space_group_type: *mut MoyoSpaceGroupType) {
    if space_group_type.is_null() {
        return;
    }
    unsafe {
        let space_group_type = Box::from_raw(space_group_type);
        free_cstring(space_group_type.hm_short);
        free_cstring(space_group_type.hm_full);
        free_cstring(space_group_type.arithmetic_symbol);
        free_cstring(space_group_type.geometric_crystal_class);
        free_cstring(space_group_type.crystal_system);
        free_cstring(space_group_type.bravais_class);
        free_cstring(space_group_type.lattice_system);
        free_cstring(space_group_type.crystal_family);
    }
}

/// Layer-group type information, created by `moyo_layer_group_type_new` and
/// freed by `moyo_layer_group_type_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoLayerGroupType {
    /// Layer-group number (1 - 80).
    pub number: i32,
    /// Hermann-Mauguin symbol in short notation.
    pub hm_short: *const c_char,
    /// Hermann-Mauguin symbol in full notation.
    pub hm_full: *const c_char,
    /// Number for layer arithmetic crystal classes (1 - 43).
    pub arithmetic_number: i32,
    /// Symbol for the layer arithmetic crystal class.
    pub arithmetic_symbol: *const c_char,
    /// Geometric crystal class. Cubic classes never occur for layer groups.
    pub geometric_crystal_class: *const c_char,
    /// Bravais class for the layer group's 2D lattice
    /// ("mp", "op", "oc", "tp", or "hp").
    pub bravais_class: *const c_char,
    /// Lattice system ("Oblique", "Rectangular", "Square", or "Hexagonal").
    pub lattice_system: *const c_char,
}

/// Look up the layer-group type information for layer-group `number` (1 - 80).
/// Returns NULL if `number` is out of range.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_layer_group_type_new(number: i32) -> *mut MoyoLayerGroupType {
    let Some(hall_number) = LayerSetting::default().hall_number(number) else {
        return std::ptr::null_mut();
    };
    let Some(lhs_entry) = layer_hall_symbol_entry(hall_number) else {
        return std::ptr::null_mut();
    };
    let Some(acc_entry) = layer_arithmetic_crystal_class_entry(lhs_entry.arithmetic_number) else {
        return std::ptr::null_mut();
    };
    let lattice_system = acc_entry.layer_lattice_system();

    Box::into_raw(Box::new(MoyoLayerGroupType {
        number,
        hm_short: leak_cstring(lhs_entry.hm_short.to_string()),
        hm_full: leak_cstring(lhs_entry.hm_full.to_string()),
        arithmetic_number: lhs_entry.arithmetic_number,
        arithmetic_symbol: leak_cstring(acc_entry.symbol.to_string()),
        geometric_crystal_class: leak_cstring(acc_entry.geometric_crystal_class.to_string()),
        bravais_class: leak_cstring(acc_entry.layer_bravais_class.to_string()),
        lattice_system: leak_cstring(lattice_system.to_string()),
    }))
}

/// Free a layer-group type created by `moyo_layer_group_type_new`.
/// Passing NULL is a no-op.
///
/// # Safety
/// `layer_group_type` must be a pointer returned by `moyo_layer_group_type_new`
/// that has not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_layer_group_type_free(layer_group_type: *mut MoyoLayerGroupType) {
    if layer_group_type.is_null() {
        return;
    }
    unsafe {
        let layer_group_type = Box::from_raw(layer_group_type);
        free_cstring(layer_group_type.hm_short);
        free_cstring(layer_group_type.hm_full);
        free_cstring(layer_group_type.arithmetic_symbol);
        free_cstring(layer_group_type.geometric_crystal_class);
        free_cstring(layer_group_type.bravais_class);
        free_cstring(layer_group_type.lattice_system);
    }
}
