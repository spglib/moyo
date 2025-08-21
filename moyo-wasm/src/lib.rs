#![cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;
use moyo::MoyoDataset;

/// Serialize errors as JS exceptions.
#[wasm_bindgen]
pub fn analyze_cell(cell_json: &str, symprec: f64, setting: &str) -> Result<JsValue, JsValue> {
    let cell: Cell = serde_json::from_str(cell_json)
        .map_err(|err| JsValue::from_str(&format!("failed to parse cell json: {}", err)))?;

    let setting = match setting {
        "Standard" | "standard" | "std" => Setting::Standard,
        "Spglib" | "spglib" => Setting::Spglib,
        _ => Setting::Standard,
    };

    let dataset = MoyoDataset::new(&cell, symprec, AngleTolerance::Default, setting)
        .map_err(|err| JsValue::from_str(&format!("moyo error: {:?}", err)))?;
    serde_wasm_bindgen::to_value(&dataset)
        .map_err(|err| JsValue::from_str(&format!("failed to serialize dataset: {}", err)))
}
