use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use crate::base::PyMoyoError;
use moyo::base::MoyoError;
use moyo::data::{LayerArithmeticNumber, layer_arithmetic_crystal_class_entry};

/// Layer arithmetic crystal class information.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "LayerArithmeticCrystalClass", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyLayerArithmeticCrystalClass {
    /// Number for layer arithmetic crystal classes (1 - 43).
    #[pyo3(get)]
    pub arithmetic_number: LayerArithmeticNumber,
    /// Symbol for the layer arithmetic crystal class (e.g. ``"p1"``, ``"c2/m11"``).
    #[pyo3(get)]
    pub symbol: &'static str,
    /// Geometric crystal class. Cubic classes never occur for layer groups.
    ///
    /// See <https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs>
    /// for string values.
    #[pyo3(get)]
    pub geometric_crystal_class: String,
    /// Bravais class for the layer group's 2D lattice (one of
    /// ``"mp"``, ``"op"``, ``"oc"``, ``"tp"``, ``"hp"``).
    #[pyo3(get)]
    pub layer_bravais_class: String,
    /// Layer lattice system (one of ``"Oblique"``, ``"Rectangular"``,
    /// ``"Square"``, ``"Hexagonal"``).
    #[pyo3(get)]
    pub layer_lattice_system: String,
}

#[pymethods]
impl PyLayerArithmeticCrystalClass {
    #[new]
    pub fn new(arithmetic_number: LayerArithmeticNumber) -> Result<Self, PyMoyoError> {
        let entry = layer_arithmetic_crystal_class_entry(arithmetic_number)
            .ok_or(MoyoError::UnknownArithmeticNumberError)?;
        let layer_lattice_system = entry.layer_lattice_system();
        Ok(Self {
            arithmetic_number: entry.arithmetic_number,
            symbol: entry.symbol,
            geometric_crystal_class: entry.geometric_crystal_class.to_string(),
            layer_bravais_class: entry.layer_bravais_class.to_string(),
            layer_lattice_system: layer_lattice_system.to_string(),
        })
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
