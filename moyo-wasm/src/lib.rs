#![cfg(target_arch = "wasm32")]
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::MoyoDataset as InternalMoyoDataset;
use moyo::base::{AngleTolerance as InternalAngleTolerance, Cell, MoyoError};
use moyo::data::{
    Centering, CrystalFamily, CrystalSystem, LatticeSystem, Setting,
    arithmetic_crystal_class_entry, hall_symbol_entry as core_hall_symbol_entry,
};
use serde::{Deserialize, Serialize};

// Generic array conversion helpers
fn to_array9<T: Copy + Into<f64>>(slice: &[T]) -> [f64; 9] {
    core::array::from_fn(|i| slice[i].into())
}

fn to_array3<T: Copy + Into<f64>>(slice: &[T]) -> [f64; 3] {
    [slice[0].into(), slice[1].into(), slice[2].into()]
}

fn parse_setting(setting: &str) -> Setting {
    if setting.to_lowercase() == "spglib" {
        Setting::Spglib
    } else {
        Setting::Standard
    }
}

fn convert_angle_tolerance(tol: InternalAngleTolerance) -> AngleTolerance {
    match tol {
        InternalAngleTolerance::Radian(x) => AngleTolerance::Radian(x),
        InternalAngleTolerance::Default => AngleTolerance::Default,
    }
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
pub struct MoyoCell {
    pub lattice: Lattice,
    pub positions: Vec<[f64; 3]>,
    pub numbers: Vec<i32>,
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
pub struct Lattice {
    pub basis: [f64; 9],
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
pub struct MoyoOperation {
    pub rotation: [f64; 9],
    pub translation: [f64; 3],
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
#[serde(tag = "type", content = "value")]
pub enum AngleTolerance {
    Radian(f64),
    Default,
}

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

fn map_moyo_err(err: MoyoError) -> JsValue {
    JsValue::from_str(&format!("moyo error: {:?}", err))
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
pub struct MoyoDataset {
    pub number: i32,
    pub hall_number: i32,
    pub hm_symbol: String,
    pub operations: Vec<MoyoOperation>,
    pub orbits: Vec<usize>,
    pub wyckoffs: Vec<String>,
    pub site_symmetry_symbols: Vec<String>,
    pub std_cell: MoyoCell,
    pub std_linear: [f64; 9],
    pub std_origin_shift: [f64; 3],
    pub std_rotation_matrix: [f64; 9],
    pub pearson_symbol: String,
    pub prim_std_cell: MoyoCell,
    pub prim_std_linear: [f64; 9],
    pub prim_std_origin_shift: [f64; 3],
    pub mapping_std_prim: Vec<usize>,
    pub symprec: f64,
    pub angle_tolerance: AngleTolerance,
}

// Inline conversions; keep file short and explicit

impl From<&moyo::base::Operation> for MoyoOperation {
    fn from(op: &moyo::base::Operation) -> Self {
        Self {
            rotation: to_array9(op.rotation.as_slice()),
            translation: to_array3(op.translation.as_slice()),
        }
    }
}

impl From<InternalMoyoDataset> for MoyoDataset {
    fn from(ds: InternalMoyoDataset) -> Self {
        let create_cell = |cell: &moyo::base::Cell| MoyoCell {
            lattice: Lattice {
                basis: to_array9(cell.lattice.basis.as_slice()),
            },
            positions: cell
                .positions
                .iter()
                .map(|p| to_array3(p.as_slice()))
                .collect(),
            numbers: cell.numbers.clone(),
        };

        Self {
            number: ds.number,
            hall_number: ds.hall_number,
            hm_symbol: ds.hm_symbol,
            operations: ds.operations.iter().map(MoyoOperation::from).collect(),
            orbits: ds.orbits,
            wyckoffs: ds.wyckoffs.iter().map(|c| c.to_string()).collect(),
            site_symmetry_symbols: ds.site_symmetry_symbols,
            std_cell: create_cell(&ds.std_cell),
            std_linear: to_array9(ds.std_linear.as_slice()),
            std_origin_shift: to_array3(ds.std_origin_shift.as_slice()),
            std_rotation_matrix: to_array9(ds.std_rotation_matrix.as_slice()),
            pearson_symbol: ds.pearson_symbol,
            prim_std_cell: create_cell(&ds.prim_std_cell),
            prim_std_linear: to_array9(ds.prim_std_linear.as_slice()),
            prim_std_origin_shift: to_array3(ds.prim_std_origin_shift.as_slice()),
            mapping_std_prim: ds.mapping_std_prim,
            symprec: ds.symprec,
            angle_tolerance: convert_angle_tolerance(ds.angle_tolerance),
        }
    }
}

/// Return a strongly-typed DTO; wasm-bindgen + tsify will emit .d.ts based on these Rust types
#[wasm_bindgen]
pub fn analyze_cell(cell_json: &str, symprec: f64, setting: &str) -> Result<MoyoDataset, JsValue> {
    // First, check if the input is valid JSON at all
    let json_value: serde_json::Value = serde_json::from_str(cell_json)
        .map_err(|err| JsValue::from_str(&format!("Input is not valid JSON: {}", err)))?;

    // Validate required fields exist
    if let Some(obj) = json_value.as_object() {
        if !obj.contains_key("lattice") {
            return Err(JsValue::from_str(
                "Missing required field 'lattice' in JSON",
            ));
        }
        if !obj.contains_key("positions") {
            return Err(JsValue::from_str(
                "Missing required field 'positions' in JSON",
            ));
        }
        if !obj.contains_key("numbers") {
            return Err(JsValue::from_str(
                "Missing required field 'numbers' in JSON",
            ));
        }

        // Check lattice structure
        if let Some(lattice) = obj.get("lattice") {
            if let Some(lattice_obj) = lattice.as_object() {
                if !lattice_obj.contains_key("basis") {
                    return Err(JsValue::from_str(
                        "Missing required field 'lattice.basis' in JSON",
                    ));
                }
            } else {
                return Err(JsValue::from_str("Field 'lattice' must be an object"));
            }
        }
    } else {
        return Err(JsValue::from_str("Input must be a JSON object"));
    }

    // Now, try to deserialize into a Cell
    let cell: Cell = serde_json::from_value(json_value).map_err(|err| {
        JsValue::from_str(&format!(
            "JSON does not match expected Cell structure: {}",
            err
        ))
    })?;

    let setting = parse_setting(setting);

    let dataset = InternalMoyoDataset::new(
        &cell,
        symprec,
        InternalAngleTolerance::default(),
        setting,
        true,
    )
    .map_err(|err| JsValue::from_str(&format!("moyo error: {:?}", err)))?;
    Ok(MoyoDataset::from(dataset))
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
