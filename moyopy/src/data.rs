use pyo3::prelude::*;
use pyo3::types::PyType;

use moyo::data::Setting;

#[derive(Debug, Clone)]
#[pyclass(name = "Setting")]
#[pyo3(module = "moyopy")]
pub struct PySetting(Setting);

#[pymethods]
impl PySetting {
    #[classmethod]
    pub fn spglib(_cls: &PyType) -> PyResult<Self> {
        Ok(Self(Setting::Spglib))
    }

    #[classmethod]
    pub fn standard(_cls: &PyType) -> PyResult<Self> {
        Ok(Self(Setting::Standard))
    }

    #[classmethod]
    pub fn hall_number(_cls: &PyType, hall_number: i32) -> PyResult<Self> {
        Ok(Self(Setting::HallNumber(hall_number)))
    }
}

impl From<PySetting> for Setting {
    fn from(setting: PySetting) -> Self {
        setting.0
    }
}
