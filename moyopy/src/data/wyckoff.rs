use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;

use moyo::data::{HallNumber, WyckoffPosition};

/// A Wyckoff position of a space group.
#[derive(Debug, Clone)]
#[pyclass(name = "WyckoffPosition", frozen, skip_from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyWyckoffPosition(WyckoffPosition);

#[derive(Serialize)]
struct PyWyckoffPositionRepr {
    hall_number: HallNumber,
    multiplicity: usize,
    letter: char,
    site_symmetry: String,
    coordinates: String,
}

impl PyWyckoffPosition {
    fn repr(&self) -> PyWyckoffPositionRepr {
        PyWyckoffPositionRepr {
            hall_number: self.0.hall_number,
            multiplicity: self.0.multiplicity,
            letter: self.0.letter,
            site_symmetry: self.0.site_symmetry.to_string(),
            coordinates: self.0.coordinates.to_string(),
        }
    }
}

#[pymethods]
impl PyWyckoffPosition {
    /// Hall symbol number of the space group this Wyckoff position belongs to.
    #[getter]
    pub fn hall_number(&self) -> HallNumber {
        self.0.hall_number
    }

    /// Site multiplicity in the conventional cell.
    #[getter]
    pub fn multiplicity(&self) -> usize {
        self.0.multiplicity
    }

    /// Wyckoff letter.
    #[getter]
    pub fn letter(&self) -> char {
        self.0.letter
    }

    /// Site-symmetry symbol.
    #[getter]
    pub fn site_symmetry(&self) -> String {
        self.0.site_symmetry.to_string()
    }

    /// Representative coordinates in the conventional setting.
    #[getter]
    pub fn coordinates(&self) -> String {
        self.0.coordinates.to_string()
    }

    fn __repr__(&self) -> String {
        self.serialize_json()
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }

    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.repr()).expect("Serialization should not fail")
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj =
                pythonize(py, &self.repr()).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }
}

impl From<WyckoffPosition> for PyWyckoffPosition {
    fn from(wyckoff: WyckoffPosition) -> Self {
        PyWyckoffPosition(wyckoff)
    }
}
