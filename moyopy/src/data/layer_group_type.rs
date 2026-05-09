use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::base::MoyoError;
use moyo::data::{
    LayerArithmeticNumber, LayerNumber, LayerSetting, layer_arithmetic_crystal_class_entry,
    layer_hall_symbol_entry,
};

use crate::base::PyMoyoError;

/// Layer-group type information.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "LayerGroupType", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyLayerGroupType {
    // Layer-group type
    /// Layer-group number (1 - 80).
    #[pyo3(get)]
    number: LayerNumber,
    /// Hermann-Mauguin symbol in short notation.
    #[pyo3(get)]
    hm_short: &'static str,
    /// Hermann-Mauguin symbol in full notation.
    #[pyo3(get)]
    hm_full: &'static str,
    // Layer arithmetic crystal class
    /// Number for layer arithmetic crystal classes (1 - 43).
    #[pyo3(get)]
    arithmetic_number: LayerArithmeticNumber,
    /// Symbol for the layer arithmetic crystal class.
    #[pyo3(get)]
    arithmetic_symbol: &'static str,
    // Other classifications
    /// Geometric crystal class. Cubic classes never occur for layer groups.
    ///
    /// See <https://github.com/spglib/moyo/blob/main/moyo/src/data/classification.rs>
    /// for string values.
    #[pyo3(get)]
    geometric_crystal_class: String,
    /// Bravais class for the layer group's 2D lattice (one of
    /// ``"mp"``, ``"op"``, ``"oc"``, ``"tp"``, ``"hp"``).
    #[pyo3(get)]
    bravais_class: String,
    /// Lattice system (one of ``"Oblique"``, ``"Rectangular"``,
    /// ``"Square"``, ``"Hexagonal"``).
    #[pyo3(get)]
    lattice_system: String,
}

#[pymethods]
impl PyLayerGroupType {
    #[new]
    pub fn new(number: LayerNumber) -> Result<Self, PyMoyoError> {
        let layer_hall_number = LayerSetting::default()
            .hall_number(number)
            .ok_or(MoyoError::UnknownNumberError)?;
        let lhs_entry =
            layer_hall_symbol_entry(layer_hall_number).ok_or(MoyoError::UnknownHallNumberError)?;

        let arithmetic_number = lhs_entry.arithmetic_number;
        let acc_entry = layer_arithmetic_crystal_class_entry(arithmetic_number)
            .ok_or(MoyoError::UnknownArithmeticNumberError)?;
        let lattice_system = acc_entry.layer_lattice_system();

        Ok(Self {
            // Layer-group type
            number,
            hm_short: lhs_entry.hm_short,
            hm_full: lhs_entry.hm_full,
            // Layer arithmetic crystal class
            arithmetic_number,
            arithmetic_symbol: acc_entry.symbol,
            // Other classifications
            geometric_crystal_class: acc_entry.geometric_crystal_class.to_string(),
            bravais_class: acc_entry.layer_bravais_class.to_string(),
            lattice_system: lattice_system.to_string(),
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
    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self).expect("Serialization should not fail")
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }
}
