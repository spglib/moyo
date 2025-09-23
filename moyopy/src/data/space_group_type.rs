use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::base::MoyoError;
use moyo::data::{
    ArithmeticNumber, CrystalFamily, CrystalSystem, LatticeSystem, Number, Setting,
    arithmetic_crystal_class_entry, hall_symbol_entry,
};

use crate::base::PyMoyoError;

#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "SpaceGroupType", frozen)]
pub struct PySpaceGroupType {
    // Space group type
    #[pyo3(get)]
    /// ITA number for space group types (1 - 230)
    number: Number,
    /// Hermann-Mauguin symbol in short notation
    #[pyo3(get)]
    hm_short: &'static str,
    /// Hermann-Mauguin symbol in full notation
    #[pyo3(get)]
    hm_full: &'static str,
    // Arithmetic crystal system
    /// Number for arithmetic crystal classes (1 - 73)
    #[pyo3(get)]
    arithmetic_number: ArithmeticNumber,
    /// Symbol for arithmetic crystal class
    #[pyo3(get)]
    arithmetic_symbol: &'static str,
    // Other classifications
    /// Geometric crystal class
    #[pyo3(get)]
    geometric_crystal_class: String,
    /// Crystal system
    #[pyo3(get)]
    crystal_system: String,
    /// Bravais class
    #[pyo3(get)]
    bravais_class: String,
    /// Lattice system
    #[pyo3(get)]
    lattice_system: String,
    /// Crystal family
    #[pyo3(get)]
    crystal_family: String,
}

#[pymethods]
impl PySpaceGroupType {
    #[new]
    pub fn new(number: Number) -> Result<Self, PyMoyoError> {
        let ita_hall_number = Setting::Standard
            .hall_number(number)
            .ok_or(MoyoError::UnknownNumberError)?;
        let ita_hall_symbol = hall_symbol_entry(ita_hall_number).unwrap();

        let arithmetic_number = ita_hall_symbol.arithmetic_number;
        let arithmetic_entry = arithmetic_crystal_class_entry(arithmetic_number).unwrap();

        let geometric_crystal_class = arithmetic_entry.geometric_crystal_class;
        let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);

        let bravais_class = arithmetic_entry.bravais_class;
        let lattice_system = LatticeSystem::from_bravais_class(bravais_class);
        let crystal_family = CrystalFamily::from_lattice_system(lattice_system);

        Ok(Self {
            // Space group type
            number,
            hm_short: ita_hall_symbol.hm_short,
            hm_full: ita_hall_symbol.hm_full,
            // Arithmetic crystal system
            arithmetic_number,
            arithmetic_symbol: arithmetic_entry.symbol,
            // Other classifications
            geometric_crystal_class: geometric_crystal_class.to_string(),
            crystal_system: crystal_system.to_string(),
            bravais_class: bravais_class.to_string(),
            lattice_system: lattice_system.to_string(),
            crystal_family: crystal_family.to_string(),
        })
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
        serde_json::to_string(&self).expect("Serialization should not fail")
    }

    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }
}
