use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{Cell, Lattice};
use moyo::utils::to_vector3;

/// A crystal structure: a lattice plus fractional positions and atomic numbers.
// Note: `PyCell` is reserved by pyo3, so we name the wrapper `PyStructure`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "Cell", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyStructure(Cell);

#[pymethods]
impl PyStructure {
    /// Create a new ``Cell``.
    ///
    /// Parameters
    /// ----------
    /// basis : list[list[float]]
    ///     Row-wise basis vectors of the lattice. ``basis[i]`` is the i-th basis vector.
    /// positions : list[list[float]]
    ///     Fractional coordinates of each site.
    /// numbers : list[int]
    ///     Atomic number of each site. Must have the same length as ``positions``.
    #[new]
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
        let positions = positions.iter().map(to_vector3).collect::<Vec<_>>();
        let cell = Cell::new(lattice, positions, numbers);

        Ok(Self(cell))
    }

    /// Row-wise basis vectors of the lattice.
    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        self.0.lattice.basis_as_array()
    }

    /// Fractional coordinates of each site.
    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0.positions_as_arrays()
    }

    /// Atomic number of each site.
    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.numbers.clone()
    }

    /// Number of atoms in the cell.
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
