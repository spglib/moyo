use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::data::{
    get_magnetic_space_group_type, ConstructType, MagneticSpaceGroupType, Number, UNINumber,
};

#[derive(Debug, Clone)]
#[pyclass(name = "MagneticSpaceGroupType", frozen)]
pub struct PyMagneticSpaceGroupType(MagneticSpaceGroupType);

#[pymethods]
impl PyMagneticSpaceGroupType {
    #[new]
    pub fn new(uni_number: UNINumber) -> Result<Self, PyErr> {
        let magnetic_space_group_type = get_magnetic_space_group_type(uni_number).ok_or(
            PyValueError::new_err(format!("Unknown uni_number: {}", uni_number)),
        )?;
        Ok(PyMagneticSpaceGroupType(magnetic_space_group_type))
    }

    #[getter]
    pub fn uni_number(&self) -> UNINumber {
        self.0.uni_number
    }

    #[getter]
    pub fn litvin_number(&self) -> i32 {
        self.0.litvin_number
    }

    #[getter]
    pub fn bns_number(&self) -> String {
        self.0.bns_number.to_string()
    }

    #[getter]
    pub fn og_number(&self) -> String {
        self.0.og_number.to_string()
    }

    #[getter]
    pub fn number(&self) -> Number {
        self.0.number
    }

    #[getter]
    pub fn construct_type(&self) -> i32 {
        match self.0.construct_type {
            ConstructType::Type1 => 1,
            ConstructType::Type2 => 2,
            ConstructType::Type3 => 3,
            ConstructType::Type4 => 4,
        }
    }
}
