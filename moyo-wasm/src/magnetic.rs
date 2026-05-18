use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::data::{
    ConstructType, get_magnetic_space_group_type,
    magnetic_hall_symbol_entry as core_magnetic_hall_symbol_entry,
};

use crate::common::{MoyoMagneticOperation, map_moyo_err};

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoMagneticHallSymbolEntry {
    pub uni_number: i32,
    pub magnetic_hall_symbol: String,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi)]
pub struct MoyoMagneticSpaceGroupType {
    pub uni_number: i32,
    pub litvin_number: i32,
    pub bns_number: String,
    pub og_number: String,
    pub number: i32,
    pub construct_type: i32,
}

fn construct_type_to_i32(ct: ConstructType) -> i32 {
    match ct {
        ConstructType::Type1 => 1,
        ConstructType::Type2 => 2,
        ConstructType::Type3 => 3,
        ConstructType::Type4 => 4,
    }
}

/// Return magnetic symmetry operations for the given `uni_number` (1 - 1651).
#[wasm_bindgen]
pub fn magnetic_operations_from_uni_number(
    uni_number: i32,
    primitive: bool,
) -> Result<Vec<MoyoMagneticOperation>, JsValue> {
    let ops = moyo::data::magnetic_operations_from_uni_number(uni_number, primitive)
        .map_err(map_moyo_err)?;
    Ok(ops.iter().map(MoyoMagneticOperation::from).collect())
}

/// Return the magnetic Hall-symbol entry for the given `uni_number` (1 - 1651).
#[wasm_bindgen]
pub fn magnetic_hall_symbol_entry(uni_number: i32) -> Result<MoyoMagneticHallSymbolEntry, JsValue> {
    let entry = core_magnetic_hall_symbol_entry(uni_number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown uni_number: {}", uni_number)))?;
    Ok(MoyoMagneticHallSymbolEntry {
        uni_number: entry.uni_number,
        magnetic_hall_symbol: entry.magnetic_hall_symbol.to_string(),
    })
}

/// Return the magnetic-space-group-type description for the given `uni_number` (1 - 1651).
#[wasm_bindgen]
pub fn magnetic_space_group_type(uni_number: i32) -> Result<MoyoMagneticSpaceGroupType, JsValue> {
    let msgt = get_magnetic_space_group_type(uni_number)
        .ok_or_else(|| JsValue::from_str(&format!("unknown uni_number: {}", uni_number)))?;
    Ok(MoyoMagneticSpaceGroupType {
        uni_number: msgt.uni_number,
        litvin_number: msgt.litvin_number,
        bns_number: msgt.bns_number.to_string(),
        og_number: msgt.og_number.to_string(),
        number: msgt.number,
        construct_type: construct_type_to_i32(msgt.construct_type),
    })
}
