#![cfg(target_arch = "wasm32")]
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::base::{AngleTolerance as InternalAngleTolerance, Cell, Position};
use moyo::data::Setting;
use moyo::MoyoDataset as InternalMoyoDataset;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

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

/// Flatten a nalgebra Matrix3<f64> into 9 f64 values in COLUMN-MAJOR order.
/// Explicit indices are used to avoid row/column confusion (nalgebra stores column-wise).
fn mat3_to_arr9_col_major(m: &Matrix3<f64>) -> [f64; 9] {
    [
        m[(0, 0)],
        m[(1, 0)],
        m[(2, 0)],
        m[(0, 1)],
        m[(1, 1)],
        m[(2, 1)],
        m[(0, 2)],
        m[(1, 2)],
        m[(2, 2)],
    ]
}

/// Flatten a nalgebra Matrix3<i32> into 9 f64 values in COLUMN-MAJOR order.
/// Explicit indices are used to avoid row/column confusion (nalgebra stores column-wise).
fn mat3i_to_arr9_col_major(m: &nalgebra::Matrix3<i32>) -> [f64; 9] {
    [
        m[(0, 0)] as f64,
        m[(1, 0)] as f64,
        m[(2, 0)] as f64,
        m[(0, 1)] as f64,
        m[(1, 1)] as f64,
        m[(2, 1)] as f64,
        m[(0, 2)] as f64,
        m[(1, 2)] as f64,
        m[(2, 2)] as f64,
    ]
}

fn vec3_to_arr3(v: &Vector3<f64>) -> [f64; 3] {
    [v[0], v[1], v[2]]
}

fn positions_to_vec3(ps: &[Position]) -> Vec<[f64; 3]> {
    ps.iter().map(|p| [p[0], p[1], p[2]]).collect()
}

fn lattice_to_dto(lattice: &moyo::base::Lattice) -> Lattice {
    Lattice {
        basis: mat3_to_arr9_col_major(&lattice.basis),
    }
}

impl From<&moyo::base::Operation> for MoyoOperation {
    fn from(op: &moyo::base::Operation) -> Self {
        MoyoOperation {
            rotation: mat3i_to_arr9_col_major(&op.rotation),
            translation: vec3_to_arr3(&op.translation),
        }
    }
}

impl From<InternalMoyoDataset> for MoyoDataset {
    fn from(ds: InternalMoyoDataset) -> Self {
        let std_cell = MoyoCell {
            lattice: lattice_to_dto(&ds.std_cell.lattice),
            positions: positions_to_vec3(&ds.std_cell.positions),
            numbers: ds.std_cell.numbers.clone(),
        };
        let prim_std_cell = MoyoCell {
            lattice: lattice_to_dto(&ds.prim_std_cell.lattice),
            positions: positions_to_vec3(&ds.prim_std_cell.positions),
            numbers: ds.prim_std_cell.numbers.clone(),
        };
        let operations: Vec<MoyoOperation> = ds
            .operations
            .iter()
            .map(|op| MoyoOperation::from(op))
            .collect();
        let angle_tolerance = match ds.angle_tolerance {
            InternalAngleTolerance::Radian(x) => AngleTolerance::Radian(x),
            InternalAngleTolerance::Default => AngleTolerance::Default,
        };
        let wyckoffs = ds.wyckoffs.iter().map(|c| c.to_string()).collect();

        MoyoDataset {
            number: ds.number,
            hall_number: ds.hall_number,
            operations,
            orbits: ds.orbits,
            wyckoffs,
            site_symmetry_symbols: ds.site_symmetry_symbols,
            std_cell,
            std_linear: mat3_to_arr9_col_major(&ds.std_linear),
            std_origin_shift: [
                ds.std_origin_shift[0],
                ds.std_origin_shift[1],
                ds.std_origin_shift[2],
            ],
            std_rotation_matrix: mat3_to_arr9_col_major(&ds.std_rotation_matrix),
            pearson_symbol: ds.pearson_symbol,
            prim_std_cell,
            prim_std_linear: mat3_to_arr9_col_major(&ds.prim_std_linear),
            prim_std_origin_shift: [
                ds.prim_std_origin_shift[0],
                ds.prim_std_origin_shift[1],
                ds.prim_std_origin_shift[2],
            ],
            mapping_std_prim: ds.mapping_std_prim,
            symprec: ds.symprec,
            angle_tolerance,
        }
    }
}

/// Return a strongly-typed DTO; wasm-bindgen + tsify will emit .d.ts based on these Rust types
#[wasm_bindgen]
pub fn analyze_cell(cell_json: &str, symprec: f64, setting: &str) -> Result<MoyoDataset, JsValue> {
    let cell: Cell = serde_json::from_str(cell_json)
        .map_err(|err| JsValue::from_str(&format!("failed to parse cell json: {}", err)))?;

    let setting = match setting {
        "Standard" | "standard" | "std" => Setting::Standard,
        "Spglib" | "spglib" => Setting::Spglib,
        _ => Setting::Standard,
    };

    let dataset =
        InternalMoyoDataset::new(&cell, symprec, InternalAngleTolerance::Default, setting)
            .map_err(|err| JsValue::from_str(&format!("moyo error: {:?}", err)))?;
    Ok(MoyoDataset::from(dataset))
}
