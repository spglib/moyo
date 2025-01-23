use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::base::MoyoError;

#[derive(Debug)]
#[pyclass(name = "MoyoError", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoError(MoyoError);

impl From<PyMoyoError> for PyErr {
    fn from(error: PyMoyoError) -> Self {
        PyValueError::new_err(error.0.to_string())
    }
}

impl From<MoyoError> for PyMoyoError {
    fn from(error: MoyoError) -> Self {
        PyMoyoError(error)
    }
}
