use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::Centering;
use moyo::utils::{to_3_slice, to_3x3_slice};

#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "Centering", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyCentering(pub Centering);

#[pymethods]
impl PyCentering {
    #[getter]
    pub fn order(&self) -> usize {
        self.0.order()
    }

    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.linear())
    }

    #[getter]
    pub fn lattice_points(&self) -> Vec<[f64; 3]> {
        self.0.lattice_points().iter().map(to_3_slice).collect()
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

impl From<PyCentering> for Centering {
    fn from(centering: PyCentering) -> Self {
        centering.0
    }
}

impl From<Centering> for PyCentering {
    fn from(centering: Centering) -> Self {
        Self(centering)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_centering() {
        let centering = PyCentering(Centering::C);
        assert_eq!(centering.order(), 2);
        let linear = centering.linear();
        assert_eq!(linear[0], [1, -1, 0]);
        assert_eq!(linear[1], [1, 1, 0]);
        assert_eq!(linear[2], [0, 0, 1]);
        assert_eq!(centering.lattice_points().len(), 2);
    }
}
