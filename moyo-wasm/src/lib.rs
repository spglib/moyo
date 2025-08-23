#![cfg(target_arch = "wasm32")]
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::base::{AngleTolerance as InternalAngleTolerance, Cell};
use moyo::data::Setting;
use moyo::MoyoDataset as InternalMoyoDataset;
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

    let dataset =
        InternalMoyoDataset::new(&cell, symprec, InternalAngleTolerance::Default, setting)
            .map_err(|err| JsValue::from_str(&format!("moyo error: {:?}", err)))?;
    Ok(MoyoDataset::from(dataset))
}
