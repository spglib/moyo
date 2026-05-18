use serde::{Deserialize, Serialize};
use tsify::Tsify;
use wasm_bindgen::prelude::*;

use moyo::base::MoyoError;
use moyo::utils::to_3_slice;

pub(crate) fn to_array9<T: Copy + Into<f64>>(slice: &[T]) -> [f64; 9] {
    core::array::from_fn(|i| slice[i].into())
}

pub(crate) fn map_moyo_err(err: MoyoError) -> JsValue {
    JsValue::from_str(&format!("moyo error: {:?}", err))
}

#[derive(Serialize, Deserialize, Tsify)]
#[tsify(into_wasm_abi, from_wasm_abi)]
pub struct MoyoOperation {
    pub rotation: [f64; 9],
    pub translation: [f64; 3],
}

impl From<&moyo::base::Operation> for MoyoOperation {
    fn from(op: &moyo::base::Operation) -> Self {
        Self {
            rotation: to_array9(op.rotation.as_slice()),
            translation: to_3_slice(&op.translation),
        }
    }
}
