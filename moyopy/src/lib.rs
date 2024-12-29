use pyo3::prelude::*;
use std::sync::OnceLock;

pub mod base;
pub mod data;

use moyo::base::AngleTolerance;
use moyo::data::Setting;
use moyo::MoyoDataset;

use crate::base::{PyMoyoError, PyOperations, PyStructure};
use crate::data::{operations_from_number, PyHallSymbolEntry, PySetting};

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

    #[getter]
    pub fn number(&self) -> i32 {
        self.0.number
    }

    #[getter]
    pub fn hall_number(&self) -> i32 {
        self.0.hall_number
    }

    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

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
}

// https://github.com/pydantic/pydantic-core/blob/main/src/lib.rs
fn moyopy_version() -> &'static str {
    static MOYOPY_VERSION: OnceLock<String> = OnceLock::new();

    MOYOPY_VERSION.get_or_init(|| {
        let version = env!("CARGO_PKG_VERSION");
        // cargo uses "1.0-alpha1" etc. while python uses "1.0.0a1", this is not full compatibility,
        // but it's good enough for now
        // see https://docs.rs/semver/1.0.9/semver/struct.Version.html#method.parse for rust spec
        // see https://peps.python.org/pep-0440/ for python spec
        // it seems the dot after "alpha/beta" e.g. "-alpha.1" is not necessary, hence why this works
        version.replace("-alpha", "a").replace("-beta", "b")
    })
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "_moyopy")]
fn moyopy(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();

    m.add("__version__", moyopy_version())?;

    // lib
    m.add_class::<PyMoyoDataset>()?;

    // base
    m.add_class::<PyStructure>()?;
    m.add_class::<PyMoyoError>()?;
    m.add_class::<PyOperations>()?;

    // data
    m.add_class::<PyHallSymbolEntry>()?;
    m.add_class::<PySetting>()?;
    m.add_wrapped(wrap_pyfunction!(operations_from_number))?;

    Ok(())
}
