use nalgebra::Vector3;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{prelude::*, IntoPyObjectExt};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{Cell, Lattice};

// Unfortunately, "PyCell" is already reversed by pyo3...
#[derive(Debug, Clone, Serialize, Deserialize)]
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

        let lattice = Lattice::from_basis(basis);
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
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::with_gil(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    #[classmethod]
    pub fn from_dict(_cls: &Bound<'_, PyType>, obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Python::with_gil(|_| {
            depythonize::<Self>(obj).map_err(|e| {
                PyErr::new::<PyValueError, _>(format!("Deserialization failed: {}", e))
            })
        })
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

#[cfg(test)]
mod tests {
    extern crate approx;

    use super::*;
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
