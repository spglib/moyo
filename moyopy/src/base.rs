use nalgebra::{OMatrix, RowVector3, Vector3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyType;
use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, Serializer};
use serde_json;

use moyo::base::{Cell, Lattice, MoyoError, Operations};

// Unfortunately, "PyCell" is already reversed by pyo3...
#[derive(Debug, Clone)]
#[pyclass(name = "Cell", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyStructure(Cell);

#[pymethods]
impl PyStructure {
    #[new]
    /// basis: row-wise basis vectors
    pub fn new(
        basis: [[f64; 3]; 3],
        positions: Vec<[f64; 3]>,
        numbers: Vec<i32>,
    ) -> PyResult<Self> {
        if positions.len() != numbers.len() {
            return Err(PyValueError::new_err(
                "positions and numbers should be the same length",
            ));
        }

        // let lattice = Lattice::new(OMatrix::from(basis));
        let lattice = Lattice::new(OMatrix::from_rows(&[
            RowVector3::from(basis[0]),
            RowVector3::from(basis[1]),
            RowVector3::from(basis[2]),
        ]));
        let positions = positions
            .iter()
            .map(|x| Vector3::new(x[0], x[1], x[2]))
            .collect::<Vec<_>>();
        let cell = Cell::new(lattice, positions, numbers);

        Ok(Self(cell))
    }

    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        *self.0.lattice.basis.as_ref()
    }

    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0.positions.iter().map(|x| [x.x, x.y, x.z]).collect()
    }

    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.numbers.clone()
    }

    #[getter]
    pub fn num_atoms(&self) -> usize {
        self.0.num_atoms()
    }

    pub fn serialize_json(&self) -> PyResult<String> {
        serde_json::to_string(self).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "Cell(basis={:?}, positions={:?}, numbers={:?})",
            self.basis(),
            self.positions(),
            self.numbers()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<PyStructure> for Cell {
    fn from(structure: PyStructure) -> Self {
        structure.0
    }
}

impl From<Cell> for PyStructure {
    fn from(cell: Cell) -> Self {
        PyStructure(cell)
    }
}

impl Serialize for PyStructure {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        Cell::from(self.clone()).serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for PyStructure {
    fn deserialize<D>(deserializer: D) -> Result<PyStructure, D::Error>
    where
        D: Deserializer<'de>,
    {
        Cell::deserialize(deserializer).map(PyStructure::from)
    }
}

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

#[cfg(test)]
mod tests {
    extern crate approx;

    use super::PyStructure;
    use approx::assert_relative_eq;
    use serde_json;

    #[test]
    fn test_serialization() {
        let structure = PyStructure::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            vec![1, 2],
        )
        .unwrap();

        let serialized = serde_json::to_string(&structure).unwrap();
        let deserialized: PyStructure = serde_json::from_str(&serialized).unwrap();

        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(structure.basis()[i][j], deserialized.basis()[i][j]);
            }
        }
        assert_eq!(structure.positions().len(), deserialized.positions().len());
        for (actual, expect) in structure.positions().iter().zip(deserialized.positions()) {
            for i in 0..3 {
                assert_relative_eq!(actual[i], expect[i]);
            }
        }
        assert_eq!(structure.numbers(), deserialized.numbers());
    }
}
