use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::LayerSetting;

/// Preference for the Hall setting of a layer group.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "LayerSetting", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyLayerSetting(pub LayerSetting);

#[pymethods]
impl PyLayerSetting {
    /// The setting with the smallest layer Hall number for each layer group
    /// (origin choice 1 for the centrosymmetric LGs 52, 62, 64).
    #[classmethod]
    pub fn spglib(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(LayerSetting::Spglib))
    }

    /// BCS / ITE standard setting per de la Flor et al., Acta Cryst. A77,
    /// 559-571 (2021), with cell choice 1 for monoclinic-oblique LGs 5 and 7,
    /// origin choice 2 for the centrosymmetric LGs 52, 62, 64, and the
    /// ``:a`` axis labelling for monoclinic-rectangular LGs 8-18.
    #[classmethod]
    pub fn standard(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(LayerSetting::Standard))
    }

    /// Specific layer Hall number from 1 to 116.
    #[classmethod]
    pub fn hall_number(_cls: &Bound<'_, PyType>, hall_number: i32) -> PyResult<Self> {
        Ok(Self(LayerSetting::HallNumber(hall_number)))
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

impl From<PyLayerSetting> for LayerSetting {
    fn from(setting: PyLayerSetting) -> Self {
        setting.0
    }
}

impl From<LayerSetting> for PyLayerSetting {
    fn from(setting: LayerSetting) -> Self {
        Self(setting)
    }
}

impl Default for PyLayerSetting {
    fn default() -> Self {
        Self(LayerSetting::default())
    }
}
