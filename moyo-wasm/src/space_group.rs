use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::data::{
    Centering, CrystalFamily, CrystalSystem, LatticeSystem, Setting,
    arithmetic_crystal_class_entry, hall_symbol_entry as core_hall_symbol_entry,
};

use crate::common::{MoyoOperation, map_moyo_err};

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(tag = "type", content = "value")]
pub enum MoyoSetting {
    Spglib,
    Standard,
    HallNumber(i32),
}

impl From<MoyoSetting> for Setting {
    fn from(s: MoyoSetting) -> Self {
        match s {
            MoyoSetting::Spglib => Setting::Spglib,
            MoyoSetting::Standard => Setting::Standard,
            MoyoSetting::HallNumber(n) => Setting::HallNumber(n),
        }
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub enum MoyoCentering {
    P,
    A,
    B,
    C,
    I,
    R,
    F,
}

impl From<Centering> for MoyoCentering {
    fn from(c: Centering) -> Self {
        match c {
            Centering::P => Self::P,
            Centering::A => Self::A,
            Centering::B => Self::B,
            Centering::C => Self::C,
            Centering::I => Self::I,
            Centering::R => Self::R,
            Centering::F => Self::F,
        }
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoHallSymbolEntry {
    pub hall_number: i32,
    pub number: i32,
    pub arithmetic_number: i32,
    pub setting: String,
    pub hall_symbol: String,
    pub hm_short: String,
    pub hm_full: String,
    pub centering: MoyoCentering,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoArithmeticCrystalClass {
    pub arithmetic_number: i32,
    pub symbol: String,
    pub geometric_crystal_class: String,
    pub bravais_class: String,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoSpaceGroupType {
    pub number: i32,
    pub hm_short: String,
    pub hm_full: String,
    pub arithmetic_number: i32,
    pub arithmetic_symbol: String,
    pub geometric_crystal_class: String,
    pub crystal_system: String,
    pub bravais_class: String,
    pub lattice_system: String,
    pub crystal_family: String,
}

/// Return symmetry operations for the given space-group ITA `number`.
#[wasm_bindgen]
pub fn operations_from_number(
    number: i32,
    setting: MoyoSetting,
    primitive: bool,
) -> Result<Vec<MoyoOperation>, JsValue> {
    let ops = moyo::data::operations_from_number(number, setting.into(), primitive)
        .map_err(map_moyo_err)?;
    Ok(ops.iter().map(MoyoOperation::from).collect())
}

/// Return the Hall-symbol entry for the given `hall_number` (1 - 530).
#[wasm_bindgen]
pub fn hall_symbol_entry(hall_number: i32) -> Result<MoyoHallSymbolEntry, JsValue> {
    let entry = core_hall_symbol_entry(hall_number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown hall_number: {}", hall_number)))?;
    Ok(MoyoHallSymbolEntry {
        hall_number: entry.hall_number,
        number: entry.number,
        arithmetic_number: entry.arithmetic_number,
        setting: entry.setting.to_string(),
        hall_symbol: entry.hall_symbol.to_string(),
        hm_short: entry.hm_short.to_string(),
        hm_full: entry.hm_full.to_string(),
        centering: entry.centering.into(),
    })
}

/// Return the arithmetic-crystal-class entry for the given `arithmetic_number` (1 - 73).
#[wasm_bindgen]
pub fn arithmetic_crystal_class(
    arithmetic_number: i32,
) -> Result<MoyoArithmeticCrystalClass, JsValue> {
    let entry = arithmetic_crystal_class_entry(arithmetic_number).ok_or_else(|| {
        JsValue::from_str(&format!("unknown arithmetic_number: {}", arithmetic_number))
    })?;
    Ok(MoyoArithmeticCrystalClass {
        arithmetic_number: entry.arithmetic_number,
        symbol: entry.symbol.to_string(),
        geometric_crystal_class: entry.geometric_crystal_class.to_string(),
        bravais_class: entry.bravais_class.to_string(),
    })
}

/// Return the space-group-type description for the given ITA `number` (1 - 230).
#[wasm_bindgen]
pub fn space_group_type(number: i32) -> Result<MoyoSpaceGroupType, JsValue> {
    let hall_number = Setting::default()
        .hall_number(number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown space-group number: {}", number)))?;
    let hs = core_hall_symbol_entry(hall_number).unwrap();
    let arith = arithmetic_crystal_class_entry(hs.arithmetic_number).unwrap();
    let crystal_system = CrystalSystem::from_geometric_crystal_class(arith.geometric_crystal_class);
    let lattice_system = LatticeSystem::from_bravais_class(arith.bravais_class);
    let crystal_family = CrystalFamily::from_lattice_system(lattice_system);
    Ok(MoyoSpaceGroupType {
        number: hs.number,
        hm_short: hs.hm_short.to_string(),
        hm_full: hs.hm_full.to_string(),
        arithmetic_number: arith.arithmetic_number,
        arithmetic_symbol: arith.symbol.to_string(),
        geometric_crystal_class: arith.geometric_crystal_class.to_string(),
        crystal_system: crystal_system.to_string(),
        bravais_class: arith.bravais_class.to_string(),
        lattice_system: lattice_system.to_string(),
        crystal_family: crystal_family.to_string(),
    })
}
