use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::Setting;

#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "Setting", frozen)]
#[pyo3(module = "moyopy")]
pub struct PySetting(pub Setting);

#[pymethods]
impl PySetting {
    #[classmethod]
    pub fn spglib(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Spglib))
    }

    #[classmethod]
    pub fn standard(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Standard))
    }

    #[classmethod]
    pub fn hall_number(_cls: &Bound<'_, PyType>, hall_number: i32) -> PyResult<Self> {
        Ok(Self(Setting::HallNumber(hall_number)))
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __repr__(&self) -> String {
        self.serialize_json()
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }

    // ------------------------------------------------------------------------
    // Serialization
    // ------------------------------------------------------------------------
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self).expect("Serialization should not fail")
    }

    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }
}

impl From<PySetting> for Setting {
    fn from(setting: PySetting) -> Self {
        setting.0
    }
}

impl From<Setting> for PySetting {
    fn from(setting: Setting) -> Self {
        Self(setting)
    }
}
