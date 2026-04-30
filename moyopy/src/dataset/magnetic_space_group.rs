use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};

use moyo::MoyoMagneticDataset;
use moyo::base::{AngleTolerance, Collinear, NonCollinear, RotationMagneticMomentAction};
use serde::{Deserialize, Serialize};

use crate::base::{
    PyCollinearMagneticCell, PyMagneticOperations, PyMoyoError, PyNonCollinearMagneticCell,
};

/// A dataset containing magnetic symmetry information of a collinear magnetic structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoCollinearMagneticDataset", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoCollinearMagneticDataset(MoyoMagneticDataset<Collinear>);

#[pymethods]
impl PyMoyoCollinearMagneticDataset {
    /// Run magnetic symmetry analysis on a collinear magnetic structure.
    ///
    /// Parameters
    /// ----------
    /// magnetic_cell : CollinearMagneticCell
    ///     Input collinear magnetic structure.
    /// symprec : float
    ///     Symmetry search tolerance in the unit of ``magnetic_cell.basis``.
    /// angle_tolerance : float | None
    ///     Symmetry search tolerance in radians.
    /// mag_symprec : float | None
    ///     Tolerance for matching magnetic moments. If ``None``, ``symprec`` is reused.
    /// is_axial : bool
    ///     If ``True``, magnetic moments transform as axial vectors. Defaults to ``False``
    ///     for collinear moments.
    /// rotate_basis : bool
    ///     Whether to rotate the basis vectors of the input cell to those of the standardized
    ///     cell.
    #[new]
    #[pyo3(signature = (magnetic_cell, *, symprec=1e-4, angle_tolerance=None, mag_symprec=None, is_axial=false, rotate_basis=true))]
    pub fn new(
        magnetic_cell: &PyCollinearMagneticCell,
        symprec: f64,
        angle_tolerance: Option<f64>,
        mag_symprec: Option<f64>,
        is_axial: bool,
        rotate_basis: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::default()
        };
        let action = if is_axial {
            RotationMagneticMomentAction::Axial
        } else {
            RotationMagneticMomentAction::Polar
        };

        let dataset = MoyoMagneticDataset::new(
            &magnetic_cell.to_owned().into(),
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
            rotate_basis,
        )?;
        Ok(PyMoyoCollinearMagneticDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------
    /// Serial number of UNI (and BNS) symbols.
    #[getter]
    pub fn uni_number(&self) -> i32 {
        self.0.uni_number
    }

    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Magnetic symmetry operations in the input cell.
    #[getter]
    pub fn magnetic_operations(&self) -> PyMagneticOperations {
        self.0.magnetic_operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// Spglib's ``crystallographic_orbits``, not ``equivalent_atoms``.
    ///
    /// The ``i`` th atom in the input magnetic cell is equivalent to the ``orbits[i]`` th
    /// atom in the **input** magnetic cell.
    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Standardized magnetic cell.
    ///
    /// Same transformation convention as :attr:`MoyoDataset.std_cell`.
    #[getter]
    pub fn std_mag_cell(&self) -> PyCollinearMagneticCell {
        self.0.std_mag_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the standardized magnetic
    /// cell.
    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        self.0.std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the standardized magnetic
    /// cell.
    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift_as_array()
    }

    /// Rigid rotation (orthogonal matrix) applied to the lattice basis.
    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        self.0.std_rotation_matrix_as_array()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Primitive standardized magnetic cell.
    #[getter]
    pub fn prim_std_mag_cell(&self) -> PyCollinearMagneticCell {
        self.0.prim_std_mag_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the primitive standardized
    /// magnetic cell.
    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        self.0.prim_std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the primitive standardized
    /// magnetic cell.
    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift_as_array()
    }

    /// Mapping of sites in the input magnetic cell to those in the primitive standardized
    /// magnetic cell.
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

    /// Actual ``mag_symprec`` used in iterative magnetic-moment matching.
    #[getter]
    pub fn mag_symprec(&self) -> f64 {
        self.0.mag_symprec
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __str__(&self) -> String {
        format!(
            "MoyoMagneticDataset(uni_number={}, magnetic_operations=<{} magnetic operations>, orbits={:?})",
            self.0.uni_number,
            self.0.magnetic_operations.len(),
            self.0.orbits
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

/// A dataset containing magnetic symmetry information of a non-collinear magnetic structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoNonCollinearMagneticDataset", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoNonCollinearMagneticDataset(MoyoMagneticDataset<NonCollinear>);

#[pymethods]
impl PyMoyoNonCollinearMagneticDataset {
    /// Run magnetic symmetry analysis on a non-collinear magnetic structure.
    ///
    /// Parameters
    /// ----------
    /// magnetic_cell : NonCollinearMagneticCell
    ///     Input non-collinear magnetic structure.
    /// symprec : float
    ///     Symmetry search tolerance in the unit of ``magnetic_cell.basis``.
    /// angle_tolerance : float | None
    ///     Symmetry search tolerance in radians.
    /// mag_symprec : float | None
    ///     Tolerance for matching magnetic moments. If ``None``, ``symprec`` is reused.
    /// is_axial : bool
    ///     If ``True``, magnetic moments transform as axial vectors. Defaults to ``True`` for
    ///     non-collinear moments.
    /// rotate_basis : bool
    ///     Whether to rotate the basis vectors of the input cell to those of the standardized
    ///     cell.
    #[new]
    #[pyo3(signature = (magnetic_cell, *, symprec=1e-4, angle_tolerance=None, mag_symprec=None, is_axial=true, rotate_basis=true))]
    pub fn new(
        magnetic_cell: &PyNonCollinearMagneticCell,
        symprec: f64,
        angle_tolerance: Option<f64>,
        mag_symprec: Option<f64>,
        is_axial: bool,
        rotate_basis: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::default()
        };
        let action = if is_axial {
            RotationMagneticMomentAction::Axial
        } else {
            RotationMagneticMomentAction::Polar
        };

        let dataset = MoyoMagneticDataset::new(
            &magnetic_cell.to_owned().into(),
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
            rotate_basis,
        )?;
        Ok(PyMoyoNonCollinearMagneticDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------
    /// Serial number of UNI (and BNS) symbols.
    #[getter]
    pub fn uni_number(&self) -> i32 {
        self.0.uni_number
    }

    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Magnetic symmetry operations in the input cell.
    #[getter]
    pub fn magnetic_operations(&self) -> PyMagneticOperations {
        self.0.magnetic_operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// Spglib's ``crystallographic_orbits``, not ``equivalent_atoms``.
    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Standardized magnetic cell.
    #[getter]
    pub fn std_mag_cell(&self) -> PyNonCollinearMagneticCell {
        self.0.std_mag_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the standardized magnetic
    /// cell.
    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        self.0.std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the standardized magnetic
    /// cell.
    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift_as_array()
    }

    /// Rigid rotation (orthogonal matrix) applied to the lattice basis.
    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        self.0.std_rotation_matrix_as_array()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------
    /// Primitive standardized magnetic cell.
    #[getter]
    pub fn prim_std_mag_cell(&self) -> PyNonCollinearMagneticCell {
        self.0.prim_std_mag_cell.clone().into()
    }

    /// Linear part of the transformation from the input cell to the primitive standardized
    /// magnetic cell.
    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        self.0.prim_std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the primitive standardized
    /// magnetic cell.
    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift_as_array()
    }

    /// Mapping of sites in the input magnetic cell to those in the primitive standardized
    /// magnetic cell.
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

    /// Actual ``mag_symprec`` used in iterative magnetic-moment matching.
    #[getter]
    pub fn mag_symprec(&self) -> f64 {
        self.0.mag_symprec
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __str__(&self) -> String {
        format!(
            "MoyoMagneticDataset(uni_number={}, magnetic_operations=<{} magnetic operations>, orbits={:?})",
            self.0.uni_number,
            self.0.magnetic_operations.len(),
            self.0.orbits
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
