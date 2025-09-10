use pyo3::prelude::*;

use crate::base::PyMoyoError;
use moyo::base::MoyoError;
use moyo::data::{arithmetic_crystal_class_entry, ArithmeticNumber};

#[derive(Debug, Clone)]
#[pyclass(name = "ArithmeticCrystalClass", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyArithmeticCrystalClass {
    // Arithmetic crystal class
    #[pyo3(get)]
    /// Number for arithmetic crystal classes (1 - 73)
    pub arithmetic_number: ArithmeticNumber,
    /// Symbol for arithmetic crystal class
    #[pyo3(get)]
    pub symbol: &'static str,
    /// Geometric crystal class
    #[pyo3(get)]
    pub geometric_crystal_class: String,
    /// Bravais class
    #[pyo3(get)]
    pub bravais_class: String,
}

#[pymethods]
impl PyArithmeticCrystalClass {
    #[new]
    pub fn new(arithmetic_number: ArithmeticNumber) -> Result<Self, PyMoyoError> {
        let entry = arithmetic_crystal_class_entry(arithmetic_number)
            .ok_or(MoyoError::UnknownArithmeticNumberError)?;
        Ok(Self {
            arithmetic_number: entry.arithmetic_number,
            symbol: entry.symbol,
            geometric_crystal_class: entry.geometric_crystal_class.to_string(),
            bravais_class: entry.bravais_class.to_string(),
        })
    }

    fn __repr__(&self) -> String {
        format!("ArithmeticCrystalClass({})", self.arithmetic_number)
    }

    fn __str__(&self) -> String {
        format!("{:?}", self)
    }
}
