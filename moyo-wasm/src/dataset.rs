use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::MoyoDataset as InternalMoyoDataset;
use moyo::base::{AngleTolerance as InternalAngleTolerance, Cell};
use moyo::data::Setting;
use moyo::utils::to_3_slice;

use crate::common::{MoyoOperation, to_array9};

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

impl From<InternalMoyoDataset> for MoyoDataset {
    fn from(ds: InternalMoyoDataset) -> Self {
        let create_cell = |cell: &moyo::base::Cell| MoyoCell {
            lattice: Lattice {
                basis: to_array9(cell.lattice.basis.as_slice()),
            },
            positions: cell.positions.iter().map(to_3_slice).collect(),
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
            std_origin_shift: to_3_slice(&ds.std_origin_shift),
            std_rotation_matrix: to_array9(ds.std_rotation_matrix.as_slice()),
            pearson_symbol: ds.pearson_symbol,
            prim_std_cell: create_cell(&ds.prim_std_cell),
            prim_std_linear: to_array9(ds.prim_std_linear.as_slice()),
            prim_std_origin_shift: to_3_slice(&ds.prim_std_origin_shift),
            mapping_std_prim: ds.mapping_std_prim,
            symprec: ds.symprec,
            angle_tolerance: convert_angle_tolerance(ds.angle_tolerance),
        }
    }
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

/// Return a strongly-typed DTO; wasm-bindgen + tsify will emit .d.ts based on these Rust types
#[wasm_bindgen]
pub fn analyze_cell(cell_json: &str, symprec: f64, setting: &str) -> Result<MoyoDataset, JsValue> {
    let json_value: serde_json::Value = serde_json::from_str(cell_json)
        .map_err(|err| JsValue::from_str(&format!("Input is not valid JSON: {}", err)))?;

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
