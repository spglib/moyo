use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{MagneticOperations, Operations};
use moyo::utils::{to_3_slice, to_3x3_slice};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "Operations", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyOperations(Operations);

#[pymethods]
impl PyOperations {
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        self.0.iter().map(|x| to_3x3_slice(&x.rotation)).collect()
    }

    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0.iter().map(|x| to_3_slice(&x.translation)).collect()
    }

    #[getter]
    pub fn num_operations(&self) -> usize {
        self.0.len()
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __len__(&self) -> usize {
        self.num_operations()
    }

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
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    #[classmethod]
    pub fn from_dict(_cls: &Bound<'_, PyType>, obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Python::attach(|_| {
            depythonize::<Self>(obj).map_err(|e| {
                PyErr::new::<PyValueError, _>(format!("Deserialization failed: {}", e))
            })
        })
    }
}

impl From<PyOperations> for Operations {
    fn from(operations: PyOperations) -> Self {
        operations.0
    }
}

impl From<Operations> for PyOperations {
    fn from(operations: Operations) -> Self {
        PyOperations(operations)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MagneticOperations", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticOperations(MagneticOperations);

#[pymethods]
impl PyMagneticOperations {
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        self.0
            .iter()
            .map(|x| to_3x3_slice(&x.operation.rotation))
            .collect()
    }

    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0
            .iter()
            .map(|x| to_3_slice(&x.operation.translation))
            .collect()
    }

    #[getter]
    pub fn time_reversals(&self) -> Vec<bool> {
        self.0.iter().map(|x| x.time_reversal).collect()
    }

    #[getter]
    pub fn num_operations(&self) -> usize {
        self.0.len()
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __len__(&self) -> usize {
        self.num_operations()
    }

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
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    #[classmethod]
    pub fn from_dict(_cls: &Bound<'_, PyType>, obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Python::attach(|_| {
            depythonize::<Self>(obj).map_err(|e| {
                PyErr::new::<PyValueError, _>(format!("Deserialization failed: {}", e))
            })
        })
    }
}

impl From<PyMagneticOperations> for MagneticOperations {
    fn from(operations: PyMagneticOperations) -> Self {
        operations.0
    }
}

impl From<MagneticOperations> for PyMagneticOperations {
    fn from(operations: MagneticOperations) -> Self {
        PyMagneticOperations(operations)
    }
}
