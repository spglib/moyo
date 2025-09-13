use pyo3::prelude::*;

use moyo::base::{MagneticOperations, Operations};
use moyo::utils::{to_3_slice, to_3x3_slice};

#[derive(Debug)]
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

    fn __len__(&self) -> usize {
        self.num_operations()
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

#[derive(Debug)]
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

    fn __len__(&self) -> usize {
        self.num_operations()
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
