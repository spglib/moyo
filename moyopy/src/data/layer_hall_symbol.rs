use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use crate::base::PyMoyoError;
use moyo::base::MoyoError;
use moyo::data::{
    LayerArithmeticNumber, LayerHallNumber, LayerHallSymbolEntry, LayerNumber,
    layer_hall_symbol_entry,
};

use super::layer_centering::PyLayerCentering;

/// An entry containing layer-group information for a specified layer
/// ``hall_number``.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "LayerHallSymbolEntry", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyLayerHallSymbolEntry(pub LayerHallSymbolEntry);

#[pymethods]
impl PyLayerHallSymbolEntry {
    #[new]
    pub fn new(hall_number: LayerHallNumber) -> Result<Self, PyMoyoError> {
        let entry =
            layer_hall_symbol_entry(hall_number).ok_or(MoyoError::UnknownHallNumberError)?;
        Ok(Self(entry.clone()))
    }

    /// Sequential number for layer-group Hall settings (1 - 116).
    #[getter]
    pub fn hall_number(&self) -> LayerHallNumber {
        self.0.hall_number
    }

    /// Layer-group number (1 - 80).
    #[getter]
    pub fn number(&self) -> LayerNumber {
        self.0.number
    }

    /// Number for layer arithmetic crystal classes (1 - 43).
    #[getter]
    pub fn arithmetic_number(&self) -> LayerArithmeticNumber {
        self.0.arithmetic_number
    }

    /// Setting code (paper Table 5 axis/origin labels: ``""``, ``"a"``,
    /// ``"b"``, ``"b-ac"``, ``"c"``, ``"c1"``, ``"c2"``, ``"c3"``, ``"1"``,
    /// ``"2"``).
    #[getter]
    pub fn setting(&self) -> &str {
        self.0.setting
    }

    /// Layer Hall symbol with lowercase ``p``/``c`` lattice prefix.
    #[getter]
    pub fn hall_symbol(&self) -> &str {
        self.0.hall_symbol
    }

    /// Hermann-Mauguin symbol in short notation.
    #[getter]
    pub fn hm_short(&self) -> &str {
        self.0.hm_short
    }

    /// Hermann-Mauguin symbol in full notation.
    #[getter]
    pub fn hm_full(&self) -> &str {
        self.0.hm_full
    }

    /// Layer centering.
    #[getter]
    pub fn centering(&self) -> PyLayerCentering {
        self.0.centering.into()
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

impl From<PyLayerHallSymbolEntry> for LayerHallSymbolEntry {
    fn from(entry: PyLayerHallSymbolEntry) -> Self {
        entry.0
    }
}

impl From<LayerHallSymbolEntry> for PyLayerHallSymbolEntry {
    fn from(entry: LayerHallSymbolEntry) -> Self {
        Self(entry)
    }
}
