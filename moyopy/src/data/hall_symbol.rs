use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use crate::base::PyMoyoError;
use moyo::base::MoyoError;
use moyo::data::{ArithmeticNumber, HallNumber, HallSymbolEntry, Number, hall_symbol_entry};

use super::centering::PyCentering;

#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "HallSymbolEntry", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyHallSymbolEntry(pub HallSymbolEntry);

#[pymethods]
impl PyHallSymbolEntry {
    #[new]
    pub fn new(hall_number: HallNumber) -> Result<Self, PyMoyoError> {
        let entry = hall_symbol_entry(hall_number).ok_or(MoyoError::UnknownHallNumberError)?;
        Ok(Self(entry))
    }

    #[getter]
    pub fn hall_number(&self) -> HallNumber {
        self.0.hall_number
    }

    #[getter]
    pub fn number(&self) -> Number {
        self.0.number
    }

    #[getter]
    pub fn arithmetic_number(&self) -> ArithmeticNumber {
        self.0.arithmetic_number
    }

    #[getter]
    pub fn setting(&self) -> &str {
        self.0.setting
    }

    #[getter]
    pub fn hall_symbol(&self) -> &str {
        self.0.hall_symbol
    }

    #[getter]
    pub fn hm_short(&self) -> &str {
        self.0.hm_short
    }

    #[getter]
    pub fn hm_full(&self) -> &str {
        self.0.hm_full
    }

    #[getter]
    pub fn centering(&self) -> PyCentering {
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

impl From<PyHallSymbolEntry> for HallSymbolEntry {
    fn from(hall_symbol_entry: PyHallSymbolEntry) -> Self {
        hall_symbol_entry.0
    }
}

impl From<HallSymbolEntry> for PyHallSymbolEntry {
    fn from(hall_symbol_entry: HallSymbolEntry) -> Self {
        Self(hall_symbol_entry)
    }
}
