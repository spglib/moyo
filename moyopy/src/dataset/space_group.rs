use pyo3::prelude::*;

use moyo::base::AngleTolerance;
use moyo::data::Setting;
use moyo::MoyoDataset;

use crate::base::{PyMoyoError, PyOperations, PyStructure};
use crate::data::PySetting;

#[derive(Debug)]
#[pyclass(name = "MoyoDataset", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoDataset(MoyoDataset);

#[pymethods]
impl PyMoyoDataset {
    #[new]
    #[pyo3(signature = (cell, *, symprec=1e-4, angle_tolerance=None, setting=None))]
    pub fn new(
        cell: &PyStructure,
        symprec: f64,
        angle_tolerance: Option<f64>,
        setting: Option<PySetting>,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::Default
        };

        let setting = if let Some(setting) = setting {
            setting.into()
        } else {
            Setting::Spglib
        };

        let dataset = MoyoDataset::new(&cell.to_owned().into(), symprec, angle_tolerance, setting)?;
        Ok(PyMoyoDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    #[getter]
    pub fn number(&self) -> i32 {
        self.0.number
    }

    #[getter]
    pub fn hall_number(&self) -> i32 {
        self.0.hall_number
    }

    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    #[getter]
    pub fn wyckoffs(&self) -> Vec<char> {
        self.0.wyckoffs.clone()
    }

    #[getter]
    pub fn site_symmetry_symbols(&self) -> Vec<String> {
        self.0.site_symmetry_symbols.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    #[getter]
    pub fn std_cell(&self) -> PyStructure {
        self.0.std_cell.clone().into()
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

    #[getter]
    pub fn pearson_symbol(&self) -> String {
        self.0.pearson_symbol.clone()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    #[getter]
    pub fn prim_std_cell(&self) -> PyStructure {
        self.0.prim_std_cell.clone().into()
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
}
