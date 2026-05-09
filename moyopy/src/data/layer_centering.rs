use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::data::LayerCentering;
use moyo::utils::{to_3_slice, to_3x3_slice};

/// Centering of a layer-group conventional cell. Only ``P`` (primitive) and
/// ``C`` (rectangular-centered) occur for layer groups.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "LayerCentering", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyLayerCentering(pub LayerCentering);

#[pymethods]
impl PyLayerCentering {
    /// Order of the centering (number of lattice points per conventional cell).
    #[getter]
    pub fn order(&self) -> usize {
        self.0.order()
    }

    /// Transformation matrix from a primitive cell to the conventional cell.
    /// The aperiodic axis ``c`` is left untouched.
    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.linear())
    }

    /// Lattice points (in fractional coordinates) of the conventional cell.
    /// The ``c`` (third) component is always zero because layer-group
    /// centerings are purely in-plane.
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

impl From<PyLayerCentering> for LayerCentering {
    fn from(centering: PyLayerCentering) -> Self {
        centering.0
    }
}

impl From<LayerCentering> for PyLayerCentering {
    fn from(centering: LayerCentering) -> Self {
        Self(centering)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_layer_centering() {
        let centering = PyLayerCentering(LayerCentering::C);
        assert_eq!(centering.order(), 2);
        let linear = centering.linear();
        assert_eq!(linear[0], [1, -1, 0]);
        assert_eq!(linear[1], [1, 1, 0]);
        assert_eq!(linear[2], [0, 0, 1]);
        assert_eq!(centering.lattice_points().len(), 2);
        for p in centering.lattice_points() {
            assert_eq!(p[2], 0.0);
        }
    }
}
