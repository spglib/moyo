use pyo3::prelude::*;

use moyo::data::Setting;

#[derive(Debug, Clone)]
#[pyclass(name = "Setting")]
#[pyo3(module = "moyo")]
pub struct PySetting(Setting);

impl From<PySetting> for Setting {
    fn from(setting: PySetting) -> Self {
        setting.0
    }
}
