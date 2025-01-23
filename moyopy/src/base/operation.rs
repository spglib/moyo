use pyo3::prelude::*;

use moyo::base::Operations;

#[derive(Debug)]
#[pyclass(name = "Operations", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyOperations(Operations);

#[pymethods]
impl PyOperations {
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0
            .iter()
            .map(|x| *x.rotation.transpose().as_ref())
            .collect()
    }

    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0.iter().map(|x| *x.translation.as_ref()).collect()
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
