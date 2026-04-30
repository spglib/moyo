use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::Setting;

/// Preference for the setting of the space group.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "Setting", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PySetting(pub Setting);

#[pymethods]
impl PySetting {
    /// The setting with the smallest Hall number.
    #[classmethod]
    pub fn spglib(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Spglib))
    }

    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral, and origin
    /// choice 2 for centrosymmetric space groups.
    #[classmethod]
    pub fn standard(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Standard))
    }

    /// Specific Hall number from 1 to 530.
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
    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self).expect("Serialization should not fail")
    }

    /// Convert this object to a dictionary.
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

impl Default for PySetting {
    fn default() -> Self {
        Self(Setting::default())
    }
}
