use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};

use moyo::MoyoDataset;
use moyo::base::AngleTolerance;
use moyo::data::Setting;

use crate::base::{PyMoyoError, PyOperations, PyStructure};
use crate::data::PySetting;

/// A dataset containing symmetry information of the input crystal structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoDataset", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoDataset(MoyoDataset);

#[pymethods]
impl PyMoyoDataset {
    /// Run symmetry analysis on a crystal structure.
    ///
    /// Parameters
    /// ----------
    /// cell : Cell
    ///     Input crystal structure.
    /// symprec : float
    ///     Symmetry search tolerance in the unit of ``cell.basis``.
    /// angle_tolerance : float | None
    ///     Symmetry search tolerance in radians.
    /// setting : Setting | None
    ///     Preference for the setting of the space group.
    /// rotate_basis : bool
    ///     Whether to rotate the basis vectors of the input cell to those of the standardized
    ///     cell.
    #[new]
    #[pyo3(signature = (cell, *, symprec=1e-4, angle_tolerance=None, setting=None, rotate_basis=true))]
    pub fn new(
        cell: &PyStructure,
        symprec: f64,
        angle_tolerance: Option<f64>,
        setting: Option<PySetting>,
        rotate_basis: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::default()
        };

        let setting = if let Some(setting) = setting {
            setting.into()
        } else {
            Setting::default()
        };

        let dataset = MoyoDataset::new(
            &cell.to_owned().into(),
            symprec,
            angle_tolerance,
            setting,
            rotate_basis,
        )?;
        Ok(PyMoyoDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    /// Space group number.
    #[getter]
    pub fn number(&self) -> i32 {
        self.0.number
    }

    /// Hall symbol number.
    #[getter]
    pub fn hall_number(&self) -> i32 {
        self.0.hall_number
    }

    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Symmetry operations in the input cell.
    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// Spglib's ``crystallographic_orbits``, not ``equivalent_atoms``.
    ///
    /// The ``i`` th atom in the input cell is equivalent to the ``orbits[i]`` th atom in the
    /// **input** cell. For example, ``orbits = [0, 0, 2, 2, 2, 2]`` means the first two atoms
    /// are equivalent and the last four atoms are equivalent to each other.
    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    /// Wyckoff letters for each site in the input cell.
    #[getter]
    pub fn wyckoffs(&self) -> Vec<char> {
        self.0.wyckoffs.clone()
    }

    /// Site symmetry symbols for each site in the input cell.
    ///
    /// The orientation of the site symmetry is with respect to the standardized cell.
    #[getter]
    pub fn site_symmetry_symbols(&self) -> Vec<String> {
        self.0.site_symmetry_symbols.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    /// Standardized cell.
    ///
    /// The input cell is related to the standardized cell by ``(std_linear, std_origin_shift)``
    /// and ``std_rotation_matrix``:
    ///
    /// ```text
    /// std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ std_linear
    /// x_std = np.linalg.inv(std_linear) @ (x_input - std_origin_shift)
    /// ```
    ///
    /// ``std_rotation_matrix`` is a rigid rotation (orthogonal matrix) applied only to the
    /// Cartesian lattice basis. It does not affect fractional coordinates.
    #[getter]
    pub fn std_cell(&self) -> PyStructure {
        self.0.std_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the standardized cell.
    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        self.0.std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the standardized cell.
    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift_as_array()
    }

    /// Rigid rotation (orthogonal matrix) applied to the lattice basis.
    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        self.0.std_rotation_matrix_as_array()
    }

    /// Pearson symbol for the standardized cell.
    #[getter]
    pub fn pearson_symbol(&self) -> String {
        self.0.pearson_symbol.clone()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    /// Primitive standardized cell.
    ///
    /// Same transformation convention as the standardized cell:
    ///
    /// ```text
    /// prim_std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ prim_std_linear
    /// x_prim_std = np.linalg.inv(prim_std_linear) @ (x_input - prim_std_origin_shift)
    /// ```
    #[getter]
    pub fn prim_std_cell(&self) -> PyStructure {
        self.0.prim_std_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the primitive standardized
    /// cell.
    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        self.0.prim_std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the primitive standardized
    /// cell.
    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift_as_array()
    }

    /// Mapping of sites in the input cell to those in the primitive standardized cell.
    ///
    /// The ``i`` th atom in the input cell is mapped to the ``mapping_std_prim[i]`` th atom in
    /// the primitive standardized cell.
    #[getter]
    pub fn mapping_std_prim(&self) -> Vec<usize> {
        self.0.mapping_std_prim.clone()
    }

    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actual ``symprec`` used in iterative symmetry search.
    #[getter]
    pub fn symprec(&self) -> f64 {
        self.0.symprec
    }

    /// Actual ``angle_tolerance`` used in iterative symmetry search.
    #[getter]
    pub fn angle_tolerance(&self) -> Option<f64> {
        if let AngleTolerance::Radian(angle_tolerance) = self.0.angle_tolerance {
            Some(angle_tolerance)
        } else {
            None
        }
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __str__(&self) -> String {
        format!(
            "MoyoDataset(number={}, hall_number={}, operations=<{} operations>, orbits={:?}, wyckoffs={:?}, site_symmetry_symbols={:?})",
            self.0.number,
            self.0.hall_number,
            self.0.operations.len(),
            self.0.orbits,
            self.0.wyckoffs,
            self.0.site_symmetry_symbols
        )
    }

    fn __repr__(&self) -> String {
        self.serialize_json()
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
