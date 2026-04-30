use pyo3::exceptions::PyValueError;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::{
    ConstructType, MagneticSpaceGroupType, Number, UNINumber, get_magnetic_space_group_type,
};

/// Magnetic space-group type information.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "MagneticSpaceGroupType", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticSpaceGroupType(MagneticSpaceGroupType);

#[pymethods]
impl PyMagneticSpaceGroupType {
    #[new]
    pub fn new(uni_number: UNINumber) -> Result<Self, PyErr> {
        let magnetic_space_group_type = get_magnetic_space_group_type(uni_number).ok_or(
            PyValueError::new_err(format!("Unknown uni_number: {}", uni_number)),
        )?;
        Ok(PyMagneticSpaceGroupType(magnetic_space_group_type))
    }

    /// Serial number of UNI (and BNS) symbols.
    #[getter]
    pub fn uni_number(&self) -> UNINumber {
        self.0.uni_number
    }

    /// Serial number in Litvin's `Magnetic group tables
    /// <https://www.iucr.org/publ/978-0-9553602-2-0>`_.
    #[getter]
    pub fn litvin_number(&self) -> i32 {
        self.0.litvin_number
    }

    /// BNS number, e.g. ``"151.32"``.
    #[getter]
    pub fn bns_number(&self) -> String {
        self.0.bns_number.to_string()
    }

    /// OG number, e.g. ``"153.4.1270"``.
    #[getter]
    pub fn og_number(&self) -> String {
        self.0.og_number.to_string()
    }

    /// ITA number for the reference space group in BNS setting.
    #[getter]
    pub fn number(&self) -> Number {
        self.0.number
    }

    /// Construct type of the magnetic space group, from 1 to 4.
    #[getter]
    pub fn construct_type(&self) -> i32 {
        match self.0.construct_type {
            ConstructType::Type1 => 1,
            ConstructType::Type2 => 2,
            ConstructType::Type3 => 3,
            ConstructType::Type4 => 4,
        }
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
