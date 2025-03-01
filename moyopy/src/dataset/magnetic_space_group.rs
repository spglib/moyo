use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{prelude::*, IntoPyObjectExt};
use pythonize::{depythonize, pythonize};

use moyo::base::{AngleTolerance, Collinear, NonCollinear, RotationMagneticMomentAction};
use moyo::MoyoMagneticDataset;
use serde::{Deserialize, Serialize};

use crate::base::{
    PyCollinearMagneticCell, PyMagneticOperations, PyMoyoError, PyNonCollinearMagneticCell,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoCollinearMagneticDataset", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoCollinearMagneticDataset(MoyoMagneticDataset<Collinear>);

#[pymethods]
impl PyMoyoCollinearMagneticDataset {
    #[new]
    #[pyo3(signature = (magnetic_cell, *, symprec=1e-4, angle_tolerance=None, mag_symprec=None, is_axial=false))]
    pub fn new(
        magnetic_cell: &PyCollinearMagneticCell,
        symprec: f64,
        angle_tolerance: Option<f64>,
        mag_symprec: Option<f64>,
        is_axial: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::Default
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
        )?;
        Ok(PyMoyoCollinearMagneticDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------

    #[getter]
    pub fn uni_number(&self) -> i32 {
        self.0.uni_number
    }

    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn magnetic_operations(&self) -> PyMagneticOperations {
        self.0.magnetic_operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------

    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn std_mag_cell(&self) -> PyCollinearMagneticCell {
        self.0.std_mag_cell.clone().into()
    }

    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.std_linear.transpose().into()
    }

    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift.into()
    }

    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.std_rotation_matrix.transpose().into()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn prim_std_mag_cell(&self) -> PyCollinearMagneticCell {
        self.0.prim_std_mag_cell.clone().into()
    }

    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.prim_std_linear.transpose().into()
    }

    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift.into()
    }

    #[getter]
    pub fn mapping_std_prim(&self) -> Vec<usize> {
        self.0.mapping_std_prim.clone()
    }

    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------

    #[getter]
    pub fn symprec(&self) -> f64 {
        self.0.symprec
    }

    #[getter]
    pub fn angle_tolerance(&self) -> Option<f64> {
        if let AngleTolerance::Radian(angle_tolerance) = self.0.angle_tolerance {
            Some(angle_tolerance)
        } else {
            None
        }
    }

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

#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoNonCollinearMagneticDataset", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoNonCollinearMagneticDataset(MoyoMagneticDataset<NonCollinear>);

#[pymethods]
impl PyMoyoNonCollinearMagneticDataset {
    #[new]
    #[pyo3(signature = (magnetic_cell, *, symprec=1e-4, angle_tolerance=None, mag_symprec=None, is_axial=true))]
    pub fn new(
        magnetic_cell: &PyNonCollinearMagneticCell,
        symprec: f64,
        angle_tolerance: Option<f64>,
        mag_symprec: Option<f64>,
        is_axial: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::Default
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
        )?;
        Ok(PyMoyoNonCollinearMagneticDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Magnetic space-group type
    // ------------------------------------------------------------------------

    #[getter]
    pub fn uni_number(&self) -> i32 {
        self.0.uni_number
    }

    // ------------------------------------------------------------------------
    // Magnetic symmetry operations in the input cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn magnetic_operations(&self) -> PyMagneticOperations {
        self.0.magnetic_operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------

    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized magnetic cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn std_mag_cell(&self) -> PyNonCollinearMagneticCell {
        self.0.std_mag_cell.clone().into()
    }

    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.std_linear.transpose().into()
    }

    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift.into()
    }

    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.std_rotation_matrix.transpose().into()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized magnetic cell
    // ------------------------------------------------------------------------

    #[getter]
    pub fn prim_std_mag_cell(&self) -> PyNonCollinearMagneticCell {
        self.0.prim_std_mag_cell.clone().into()
    }

    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        // Since nalgebra stores matrices in column-major order, we need to transpose them
        self.0.prim_std_linear.transpose().into()
    }

    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift.into()
    }

    #[getter]
    pub fn mapping_std_prim(&self) -> Vec<usize> {
        self.0.mapping_std_prim.clone()
    }

    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------

    #[getter]
    pub fn symprec(&self) -> f64 {
        self.0.symprec
    }

    #[getter]
    pub fn angle_tolerance(&self) -> Option<f64> {
        if let AngleTolerance::Radian(angle_tolerance) = self.0.angle_tolerance {
            Some(angle_tolerance)
        } else {
            None
        }
    }

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
