use pyo3::prelude::*;
use pyo3::types::PyType;

use super::base::{PyMoyoError, PyOperations};
use moyo::data::{hall_symbol_entry, HallSymbol};
use moyo::{MoyoError, Operations, Setting};

#[derive(Debug, Clone)]
#[pyclass(name = "Setting")]
#[pyo3(module = "moyopy")]
pub struct PySetting(Setting);

#[pymethods]
impl PySetting {
    #[classmethod]
    pub fn spglib(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Spglib))
    }

    #[classmethod]
    pub fn standard(_cls: &Bound<'_, PyType>) -> PyResult<Self> {
        Ok(Self(Setting::Standard))
    }

    #[classmethod]
    pub fn hall_number(_cls: &Bound<'_, PyType>, hall_number: i32) -> PyResult<Self> {
        Ok(Self(Setting::HallNumber(hall_number)))
    }
}

impl From<PySetting> for Setting {
    fn from(setting: PySetting) -> Self {
        setting.0
    }
}

#[pyfunction]
pub fn operations_from_number(
    number: i32,
    setting: Option<PySetting>,
) -> Result<PyOperations, PyMoyoError> {
    let setting = if let Some(setting) = setting {
        setting
    } else {
        PySetting(Setting::Spglib)
    };
    let hall_number = match setting.0 {
        Setting::HallNumber(hall_number) => hall_number,
        Setting::Spglib | Setting::Standard => setting.0.hall_numbers()[(number - 1) as usize],
    };
    let entry = hall_symbol_entry(hall_number);
    let hs = HallSymbol::new(entry.hall_symbol).ok_or(MoyoError::HallSymbolParsingError)?;

    let mut rotations = vec![];
    let mut translations = vec![];

    let coset = hs.traverse();
    let lattice_points = hs.centering.lattice_points();
    for t1 in lattice_points.iter() {
        for (r2, t2) in coset.rotations.iter().zip(coset.translations.iter()) {
            // (E, t1) (r2, t2) = (r2, t1 + t2)
            rotations.push(*r2);
            let t12 = (t1 + t2).map(|e| e % 1.);
            translations.push(t12);
        }
    }
    Ok(PyOperations::from(Operations::new(rotations, translations)))
}
