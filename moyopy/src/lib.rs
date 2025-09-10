use pyo3::prelude::*;
use std::sync::OnceLock;

mod utils;

pub mod base;
pub mod data;
pub mod dataset;
pub mod identify;

use crate::base::{
    PyCollinearMagneticCell, PyMagneticOperations, PyMoyoError, PyNonCollinearMagneticCell,
    PyOperations, PyStructure,
};
use crate::data::{
    operations_from_number, PyArithmeticCrystalClass, PyCentering, PyHallSymbolEntry,
    PyMagneticSpaceGroupType, PySetting, PySpaceGroupType,
};
use crate::dataset::{
    PyMoyoCollinearMagneticDataset, PyMoyoDataset, PyMoyoNonCollinearMagneticDataset,
};
use crate::identify::{PyPointGroup, PySpaceGroup};

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

    // lib
    m.add("__version__", moyopy_version())?;

    // dataset
    m.add_class::<PyMoyoDataset>()?;
    m.add_class::<PyMoyoCollinearMagneticDataset>()?;
    m.add_class::<PyMoyoNonCollinearMagneticDataset>()?;

    // base
    m.add_class::<PyStructure>()?;
    m.add_class::<PyCollinearMagneticCell>()?;
    m.add_class::<PyNonCollinearMagneticCell>()?;
    m.add_class::<PyMoyoError>()?;
    m.add_class::<PyOperations>()?;
    m.add_class::<PyMagneticOperations>()?;

    // data: Hall symbol data
    m.add_class::<PySetting>()?;
    m.add_class::<PyCentering>()?;
    m.add_class::<PyHallSymbolEntry>()?;
    // data: Group data
    m.add_class::<PySpaceGroupType>()?;
    m.add_class::<PyMagneticSpaceGroupType>()?;
    m.add_class::<PyArithmeticCrystalClass>()?;
    // data: Misc
    m.add_wrapped(wrap_pyfunction!(operations_from_number))?;

    // identify
    m.add_class::<PyPointGroup>()?;
    m.add_class::<PySpaceGroup>()?;

    Ok(())
}
