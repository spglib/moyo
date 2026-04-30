use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};
use serde_json;

use moyo::base::{Collinear, Lattice, MagneticCell, NonCollinear};
use moyo::utils::to_vector3;

#[derive(Debug, Clone, Serialize, Deserialize)]
enum MagneticCellEnum {
    CollinearMagneticCell(MagneticCell<Collinear>),
    NonCollinearMagneticCell(MagneticCell<NonCollinear>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MagneticCell", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticCell {
    magnetic_cell: MagneticCellEnum,
}

/// A crystal structure with collinear magnetic moments.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "CollinearMagneticCell", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyCollinearMagneticCell(MagneticCell<Collinear>);

#[pymethods]
impl PyCollinearMagneticCell {
    /// Create a new ``CollinearMagneticCell``.
    ///
    /// Parameters
    /// ----------
    /// basis : list[list[float]]
    ///     Row-wise basis vectors of the lattice.
    /// positions : list[list[float]]
    ///     Fractional coordinates of each site.
    /// numbers : list[int]
    ///     Atomic number of each site.
    /// magnetic_moments : list[float]
    ///     Scalar magnetic moment of each site (collinear).
    #[new]
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
        let positions = positions.iter().map(to_vector3).collect::<Vec<_>>();
        let magnetic_moments = magnetic_moments.iter().map(|m| Collinear(*m)).collect();
        let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magnetic_moments);

        Ok(Self(magnetic_cell))
    }

    /// Row-wise basis vectors of the lattice.
    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        self.0.cell.lattice.basis_as_array()
    }

    /// Fractional coordinates of each site.
    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0.cell.positions_as_arrays()
    }

    /// Atomic number of each site.
    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.cell.numbers.clone()
    }

    /// Scalar magnetic moment of each site.
    #[getter]
    pub fn magnetic_moments(&self) -> Vec<f64> {
        self.0.magnetic_moments.iter().map(|m| m.0).collect()
    }

    /// Number of atoms in the magnetic cell.
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

/// A crystal structure with non-collinear (vector) magnetic moments.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "NonCollinearMagneticCell", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyNonCollinearMagneticCell(MagneticCell<NonCollinear>);

#[pymethods]
impl PyNonCollinearMagneticCell {
    /// Create a new ``NonCollinearMagneticCell``.
    ///
    /// Parameters
    /// ----------
    /// basis : list[list[float]]
    ///     Row-wise basis vectors of the lattice.
    /// positions : list[list[float]]
    ///     Fractional coordinates of each site.
    /// numbers : list[int]
    ///     Atomic number of each site.
    /// magnetic_moments : list[list[float]]
    ///     Three-component magnetic moment vector of each site.
    #[new]
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
        let positions = positions.iter().map(to_vector3).collect::<Vec<_>>();
        let magnetic_moments = magnetic_moments
            .iter()
            .map(|m| NonCollinear(to_vector3(m)))
            .collect();
        let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magnetic_moments);

        Ok(Self(magnetic_cell))
    }

    /// Row-wise basis vectors of the lattice.
    #[getter]
    pub fn basis(&self) -> [[f64; 3]; 3] {
        self.0.cell.lattice.basis_as_array()
    }

    /// Fractional coordinates of each site.
    #[getter]
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.0.cell.positions_as_arrays()
    }

    /// Atomic number of each site.
    #[getter]
    pub fn numbers(&self) -> Vec<i32> {
        self.0.cell.numbers.clone()
    }

    /// Three-component magnetic moment vector of each site.
    #[getter]
    pub fn magnetic_moments(&self) -> Vec<[f64; 3]> {
        self.0
            .magnetic_moments
            .iter()
            .map(|m| m.as_array())
            .collect()
    }

    /// Number of atoms in the magnetic cell.
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
