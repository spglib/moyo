use nalgebra::Vector3;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{prelude::*, IntoPyObjectExt};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{Collinear, Lattice, MagneticCell, NonCollinear};

#[derive(Debug, Clone, Serialize, Deserialize)]
enum MagneticCellEnum {
    CollinearMagneticCell(MagneticCell<Collinear>),
    NonCollinearMagneticCell(MagneticCell<NonCollinear>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MagneticCell", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticCell {
    magnetic_cell: MagneticCellEnum,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "CollinearMagneticCell", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyCollinearMagneticCell(MagneticCell<Collinear>);

#[pymethods]
impl PyCollinearMagneticCell {
    #[new]
    /// basis: row-wise basis vectors
    pub fn new(
        basis: [[f64; 3]; 3],
        positions: Vec<[f64; 3]>,
        numbers: Vec<i32>,
        magnetic_moments: Vec<f64>,
    ) -> PyResult<Self> {
        if numbers.len() != positions.len() {
            return Err(PyValueError::new_err(
                "positions and numbers should be the same length",
            ));
        }
        if magnetic_moments.len() != positions.len() {
            return Err(PyValueError::new_err(
                "positions and magnetic_moments should be the same length",
            ));
        }

        let lattice = Lattice::from_basis(basis);
        let positions = positions
            .iter()
            .map(|x| Vector3::new(x[0], x[1], x[2]))
            .collect::<Vec<_>>();
        let magnetic_moments = magnetic_moments.iter().map(|m| Collinear(*m)).collect();
        let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magnetic_moments);

        Ok(Self(magnetic_cell))
    }

    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        *self.0.cell.lattice.basis.as_ref()
    }

    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0
            .cell
            .positions
            .iter()
            .map(|x| [x.x, x.y, x.z])
            .collect()
    }

    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.cell.numbers.clone()
    }

    #[getter]
    pub fn magnetic_moments(&self) -> Vec<f64> {
        self.0.magnetic_moments.iter().map(|m| m.0).collect()
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

impl From<PyCollinearMagneticCell> for MagneticCell<Collinear> {
    fn from(cell: PyCollinearMagneticCell) -> Self {
        cell.0
    }
}

impl From<MagneticCell<Collinear>> for PyCollinearMagneticCell {
    fn from(cell: MagneticCell<Collinear>) -> Self {
        PyCollinearMagneticCell(cell)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "NonCollinearMagneticCell", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyNonCollinearMagneticCell(MagneticCell<NonCollinear>);

#[pymethods]
impl PyNonCollinearMagneticCell {
    #[new]
    /// basis: row-wise basis vectors
    pub fn new(
        basis: [[f64; 3]; 3],
        positions: Vec<[f64; 3]>,
        numbers: Vec<i32>,
        magnetic_moments: Vec<[f64; 3]>,
    ) -> PyResult<Self> {
        if numbers.len() != positions.len() {
            return Err(PyValueError::new_err(
                "positions and numbers should be the same length",
            ));
        }
        if magnetic_moments.len() != positions.len() {
            return Err(PyValueError::new_err(
                "positions and magnetic_moments should be the same length",
            ));
        }

        let lattice = Lattice::from_basis(basis);
        let positions = positions
            .iter()
            .map(|x| Vector3::new(x[0], x[1], x[2]))
            .collect::<Vec<_>>();
        let magnetic_moments = magnetic_moments
            .iter()
            .map(|m| NonCollinear(Vector3::new(m[0], m[1], m[2])))
            .collect();
        let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magnetic_moments);

        Ok(Self(magnetic_cell))
    }

    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        *self.0.cell.lattice.basis.as_ref()
    }

    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0
            .cell
            .positions
            .iter()
            .map(|x| [x.x, x.y, x.z])
            .collect()
    }

    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.cell.numbers.clone()
    }

    #[getter]
    pub fn magnetic_moments(&self) -> Vec<[f64; 3]> {
        self.0
            .magnetic_moments
            .iter()
            .map(|m| [m.0[0], m.0[1], m.0[2]])
            .collect()
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

impl From<PyNonCollinearMagneticCell> for MagneticCell<NonCollinear> {
    fn from(cell: PyNonCollinearMagneticCell) -> Self {
        cell.0
    }
}

impl From<MagneticCell<NonCollinear>> for PyNonCollinearMagneticCell {
    fn from(cell: MagneticCell<NonCollinear>) -> Self {
        PyNonCollinearMagneticCell(cell)
    }
}
