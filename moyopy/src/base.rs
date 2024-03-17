use nalgebra::{OMatrix, RowVector3, Vector3};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::base::{Cell, Lattice, MoyoError, Operations};

// Unfortunately, "PyCell" is already reversed by pyo3...
#[derive(Debug, Clone)]
#[pyclass(name = "Cell")]
#[pyo3(module = "moyo")]
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

#[derive(Debug)]
#[pyclass(name = "MoyoError")]
#[pyo3(module = "moyo")]
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
#[pyclass(name = "Operations")]
#[pyo3(module = "moyo")]
pub struct PyOperations(Operations);

#[pymethods]
impl PyOperations {
    #[getter]
    pub fn rotations(&self) -> Vec<[[i32; 3]; 3]> {
        self.0.rotations.iter().map(|x| *x.as_ref()).collect()
    }

    #[getter]
    pub fn translations(&self) -> Vec<[f64; 3]> {
        self.0.translations.iter().map(|x| *x.as_ref()).collect()
    }

    #[getter]
    pub fn num_operations(&self) -> usize {
        self.0.num_operations()
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
