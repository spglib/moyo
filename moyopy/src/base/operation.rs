use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{MagneticOperations, Operations, UnimodularTransformation};

/// A list of crystallographic symmetry operations (rotation + translation).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "Operations", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyOperations(Operations);

#[pymethods]
impl PyOperations {
    /// Rotation parts of the symmetry operations.
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        self.0.iter().map(|x| x.rotation_as_array()).collect()
    }

    /// Translation parts of the symmetry operations (fractional coordinates).
    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0.iter().map(|x| x.translation_as_array()).collect()
    }

    /// Number of symmetry operations.
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
    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    /// Deserialize an object from a JSON string.
    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    /// Create an object from a dictionary.
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

/// A list of magnetic symmetry operations (rotation + translation + time reversal).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MagneticOperations", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticOperations(MagneticOperations);

#[pymethods]
impl PyMagneticOperations {
    /// Rotation parts of the magnetic symmetry operations.
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        self.0
            .iter()
            .map(|x| x.operation.rotation_as_array())
            .collect()
    }

    /// Translation parts of the magnetic symmetry operations (fractional coordinates).
    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0
            .iter()
            .map(|x| x.operation.translation_as_array())
            .collect()
    }

    /// Time-reversal flag for each magnetic symmetry operation.
    #[getter]
    pub fn time_reversals(&self) -> Vec<bool> {
        self.0.iter().map(|x| x.time_reversal).collect()
    }

    /// Number of magnetic symmetry operations.
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
    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    /// Deserialize an object from a JSON string.
    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    /// Create an object from a dictionary.
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

/// A unimodular transformation: an integer linear part with determinant +-1 plus an origin
/// shift.
#[derive(Debug, Clone)]
#[pyclass(name = "UnimodularTransformation", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyUnimodularTransformation(UnimodularTransformation);

#[derive(Serialize)]
struct PyUnimodularTransformationRepr {
    linear: [[i32; 3]; 3],
    origin_shift: [f64; 3],
}

impl PyUnimodularTransformation {
    fn repr(&self) -> PyUnimodularTransformationRepr {
        PyUnimodularTransformationRepr {
            linear: self.0.linear_as_array(),
            origin_shift: self.0.origin_shift_as_array(),
        }
    }
}

#[pymethods]
impl PyUnimodularTransformation {
    /// Linear part (integer 3x3 matrix with determinant +-1).
    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        self.0.linear_as_array()
    }

    /// Origin shift (fractional coordinates).
    #[getter]
    pub fn origin_shift(&self) -> [f64; 3] {
        self.0.origin_shift_as_array()
    }

    fn __repr__(&self) -> String {
        self.serialize_json()
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }

    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.repr()).expect("Serialization should not fail")
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj =
                pythonize(py, &self.repr()).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }
}

impl From<PyUnimodularTransformation> for UnimodularTransformation {
    fn from(transformation: PyUnimodularTransformation) -> Self {
        transformation.0
    }
}

impl From<UnimodularTransformation> for PyUnimodularTransformation {
    fn from(transformation: UnimodularTransformation) -> Self {
        PyUnimodularTransformation(transformation)
    }
}
