use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::data::{
    CrystalSystem, LayerCentering, LayerHallSymbolEntry, LayerSetting,
    iter_layer_hall_symbol_entry, iter_layer_wyckoff_positions_from_hall_number,
    layer_arithmetic_crystal_class_entry, layer_hall_symbol_entry as core_layer_hall_symbol_entry,
};

use crate::common::{MoyoOperation, map_moyo_err};

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(tag = "type", content = "value")]
pub enum MoyoLayerSetting {
    Spglib,
    Standard,
    HallNumber(i32),
}

impl From<MoyoLayerSetting> for LayerSetting {
    fn from(s: MoyoLayerSetting) -> Self {
        match s {
            MoyoLayerSetting::Spglib => LayerSetting::Spglib,
            MoyoLayerSetting::Standard => LayerSetting::Standard,
            MoyoLayerSetting::HallNumber(n) => LayerSetting::HallNumber(n),
        }
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub enum MoyoLayerCentering {
    P,
    C,
}

impl From<LayerCentering> for MoyoLayerCentering {
    fn from(c: LayerCentering) -> Self {
        match c {
            LayerCentering::P => Self::P,
            LayerCentering::C => Self::C,
        }
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoLayerHallSymbolEntry {
    pub hall_number: i32,
    pub number: i32,
    pub arithmetic_number: i32,
    pub setting: String,
    pub hall_symbol: String,
    pub hm_short: String,
    pub hm_full: String,
    pub centering: MoyoLayerCentering,
}

impl From<&LayerHallSymbolEntry> for MoyoLayerHallSymbolEntry {
    fn from(entry: &LayerHallSymbolEntry) -> Self {
        Self {
            hall_number: entry.hall_number,
            number: entry.number,
            arithmetic_number: entry.arithmetic_number,
            setting: entry.setting.to_string(),
            hall_symbol: entry.hall_symbol.to_string(),
            hm_short: entry.hm_short.to_string(),
            hm_full: entry.hm_full.to_string(),
            centering: entry.centering.into(),
        }
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoLayerWyckoffPosition {
    /// Multiplicity in the conventional layer cell.
    pub multiplicity: usize,
    /// Wyckoff letter.
    pub letter: String,
    /// Site-symmetry symbol.
    pub site_symmetry: String,
    /// Representative coordinate triplet (e.g. `x,y,z`).
    pub coordinates: String,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoLayerArithmeticCrystalClass {
    pub arithmetic_number: i32,
    pub symbol: String,
    pub geometric_crystal_class: String,
    pub bravais_class: String,
    pub lattice_system: String,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoLayerGroupType {
    pub number: i32,
    pub hall_number: i32,
    pub hm_short: String,
    pub hm_full: String,
    pub arithmetic_number: i32,
    pub arithmetic_symbol: String,
    pub geometric_crystal_class: String,
    pub crystal_system: String,
    pub bravais_class: String,
    pub lattice_system: String,
}

/// Return symmetry operations for the given layer-group `number` (1 - 80).
#[wasm_bindgen]
pub fn operations_from_layer_number(
    number: i32,
    setting: MoyoLayerSetting,
    primitive: bool,
) -> Result<Vec<MoyoOperation>, JsValue> {
    let ops = moyo::data::operations_from_layer_number(number, setting.into(), primitive)
        .map_err(map_moyo_err)?;
    Ok(ops.iter().map(MoyoOperation::from).collect())
}

/// Return the layer Hall-symbol entry for the given `hall_number` (1 - 116).
#[wasm_bindgen]
pub fn layer_hall_symbol_entry(hall_number: i32) -> Result<MoyoLayerHallSymbolEntry, JsValue> {
    let entry = core_layer_hall_symbol_entry(hall_number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown layer hall_number: {}", hall_number)))?;
    Ok(entry.into())
}

/// Return all layer Hall-symbol entries for the given layer-group `number` (1 - 80),
/// ordered by Hall number. Each entry corresponds to one setting of the layer group.
#[wasm_bindgen]
pub fn layer_hall_symbol_entries_from_number(number: i32) -> Vec<MoyoLayerHallSymbolEntry> {
    iter_layer_hall_symbol_entry()
        .filter(|entry| entry.number == number)
        .map(MoyoLayerHallSymbolEntry::from)
        .collect()
}

/// Return all Wyckoff positions for the given layer Hall number (1 - 116), ordered
/// general-position-first then descending multiplicity. Returns an empty array for
/// unknown Hall numbers.
#[wasm_bindgen]
pub fn layer_wyckoff_positions(hall_number: i32) -> Vec<MoyoLayerWyckoffPosition> {
    iter_layer_wyckoff_positions_from_hall_number(hall_number)
        .map(|wp| MoyoLayerWyckoffPosition {
            multiplicity: wp.multiplicity,
            letter: wp.letter.to_string(),
            site_symmetry: wp.site_symmetry.to_string(),
            coordinates: wp.coordinates.to_string(),
        })
        .collect()
}

/// Return the layer arithmetic-crystal-class entry for the given `arithmetic_number` (1 - 43).
#[wasm_bindgen]
pub fn layer_arithmetic_crystal_class(
    arithmetic_number: i32,
) -> Result<MoyoLayerArithmeticCrystalClass, JsValue> {
    let entry = layer_arithmetic_crystal_class_entry(arithmetic_number).ok_or_else(|| {
        JsValue::from_str(&format!(
            "unknown layer arithmetic_number: {}",
            arithmetic_number
        ))
    })?;
    let lattice_system = entry.layer_lattice_system();
    Ok(MoyoLayerArithmeticCrystalClass {
        arithmetic_number: entry.arithmetic_number,
        symbol: entry.symbol.to_string(),
        geometric_crystal_class: entry.geometric_crystal_class.to_string(),
        bravais_class: entry.layer_bravais_class.to_string(),
        lattice_system: lattice_system.to_string(),
    })
}

/// Return the layer-group-type description for the given `number` (1 - 80).
#[wasm_bindgen]
pub fn layer_group_type(number: i32) -> Result<MoyoLayerGroupType, JsValue> {
    let layer_hall_number = LayerSetting::default()
        .hall_number(number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown layer-group number: {}", number)))?;
    let lhs = core_layer_hall_symbol_entry(layer_hall_number).unwrap();
    let acc = layer_arithmetic_crystal_class_entry(lhs.arithmetic_number).unwrap();
    let crystal_system = CrystalSystem::from_geometric_crystal_class(acc.geometric_crystal_class);
    let lattice_system = acc.layer_lattice_system();
    Ok(MoyoLayerGroupType {
        number: lhs.number,
        hall_number: lhs.hall_number,
        hm_short: lhs.hm_short.to_string(),
        hm_full: lhs.hm_full.to_string(),
        arithmetic_number: acc.arithmetic_number,
        arithmetic_symbol: acc.symbol.to_string(),
        geometric_crystal_class: acc.geometric_crystal_class.to_string(),
        crystal_system: crystal_system.to_string(),
        bravais_class: acc.layer_bravais_class.to_string(),
        lattice_system: lattice_system.to_string(),
    })
}
