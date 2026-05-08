use pyo3::exceptions::PyValueError;
use pyo3::types::PyType;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::{depythonize, pythonize};
use serde::{Deserialize, Serialize};

use moyo::MoyoLayerDataset;
use moyo::base::{AngleTolerance, Cell, Lattice, LayerCell};
use moyo::data::LayerSetting;

use crate::base::{PyMoyoError, PyOperations, PyStructure};
use crate::data::PyLayerSetting;

/// Re-pack a [`LayerCell`] as a bulk [`Cell`] using its public getters,
/// so the moyopy boundary can hand the standardized layer cell back to
/// Python as a [`PyStructure`]. The bulk-vs-layer split is enforced inside
/// `moyo`; this conversion is a one-shot output bridge.
fn layer_cell_to_cell(layer: &LayerCell) -> Cell {
    Cell::new(
        Lattice {
            basis: *layer.lattice().basis(),
        },
        layer.positions().to_vec(),
        layer.numbers().to_vec(),
    )
}

/// A dataset containing layer-group symmetry information of the input crystal
/// structure (a 2D-periodic system whose third basis vector is the aperiodic
/// stacking direction).
#[derive(Debug, Clone, Serialize, Deserialize)]
#[pyclass(name = "MoyoLayerDataset", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMoyoLayerDataset(MoyoLayerDataset);

#[pymethods]
impl PyMoyoLayerDataset {
    /// Run layer-group symmetry analysis on a crystal structure.
    ///
    /// Parameters
    /// ----------
    /// cell : Cell
    ///     Input crystal structure. The third basis vector ``c`` must be the
    ///     aperiodic stacking direction and perpendicular to ``a, b``;
    ///     ``a, b`` must lie in the xy-plane. Inputs that violate this contract
    ///     are rejected with a descriptive error.
    /// symprec : float
    ///     Symmetry search tolerance in the unit of ``cell.basis``.
    /// angle_tolerance : float | None
    ///     Symmetry search tolerance in radians.
    /// setting : LayerSetting | None
    ///     Preference for the Hall setting of the layer group.
    /// rotate_basis : bool
    ///     Whether to rotate the basis vectors of the input cell to those of
    ///     the standardized cell.
    #[new]
    #[pyo3(signature = (cell, *, symprec=1e-4, angle_tolerance=None, setting=None, rotate_basis=true))]
    pub fn new(
        cell: &PyStructure,
        symprec: f64,
        angle_tolerance: Option<f64>,
        setting: Option<PyLayerSetting>,
        rotate_basis: bool,
    ) -> Result<Self, PyMoyoError> {
        let angle_tolerance = if let Some(angle_tolerance) = angle_tolerance {
            AngleTolerance::Radian(angle_tolerance)
        } else {
            AngleTolerance::default()
        };

        let setting = if let Some(setting) = setting {
            setting.into()
        } else {
            LayerSetting::default()
        };

        let dataset = MoyoLayerDataset::new(
            &cell.to_owned().into(),
            symprec,
            angle_tolerance,
            setting,
            rotate_basis,
        )?;
        Ok(PyMoyoLayerDataset(dataset))
    }

    // ------------------------------------------------------------------------
    // Identification
    // ------------------------------------------------------------------------
    /// Layer group number (1 - 80).
    #[getter]
    pub fn number(&self) -> i32 {
        self.0.number
    }

    /// Layer Hall symbol number (1 - 116).
    #[getter]
    pub fn hall_number(&self) -> i32 {
        self.0.hall_number
    }

    // ------------------------------------------------------------------------
    // Symmetry operations in the input cell
    // ------------------------------------------------------------------------
    /// Layer-group operations in the input cell.
    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

    // ------------------------------------------------------------------------
    // Site symmetry
    // ------------------------------------------------------------------------
    /// The ``i`` th atom in the input cell is equivalent to the ``orbits[i]``
    /// th atom in the input cell.
    #[getter]
    pub fn orbits(&self) -> Vec<usize> {
        self.0.orbits.clone()
    }

    /// Wyckoff letters for each site in the input cell.
    #[getter]
    pub fn wyckoffs(&self) -> Vec<char> {
        self.0.wyckoffs.clone()
    }

    /// Site symmetry symbols for each site in the input cell, oriented w.r.t.
    /// the standardized cell.
    #[getter]
    pub fn site_symmetry_symbols(&self) -> Vec<String> {
        self.0.site_symmetry_symbols.clone()
    }

    // ------------------------------------------------------------------------
    // Standardized layer cell
    // ------------------------------------------------------------------------
    /// Conventional standardized layer cell.
    ///
    /// The input cell is related to the standardized cell by
    /// ``(std_linear, std_origin_shift)`` and ``std_rotation_matrix``:
    ///
    /// ```text
    /// std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ std_linear
    /// x_std = np.linalg.inv(std_linear) @ (x_input - std_origin_shift)
    /// ```
    #[getter]
    pub fn std_cell(&self) -> PyStructure {
        layer_cell_to_cell(&self.0.std_cell).into()
    }

    /// Linear part of the transformation from the input cell to the
    /// standardized layer cell.
    #[getter]
    pub fn std_linear(&self) -> [[f64; 3]; 3] {
        self.0.std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the
    /// standardized layer cell.
    #[getter]
    pub fn std_origin_shift(&self) -> [f64; 3] {
        self.0.std_origin_shift_as_array()
    }

    /// Rigid rotation (orthogonal matrix) applied to the lattice basis when
    /// ``rotate_basis = True``, identity otherwise.
    #[getter]
    pub fn std_rotation_matrix(&self) -> [[f64; 3]; 3] {
        self.0.std_rotation_matrix_as_array()
    }

    /// Pearson symbol for the standardized layer cell. The first two
    /// characters are the 2D Bravais type (``mp``, ``op``, ``oc``, ``tp``,
    /// or ``hp``).
    #[getter]
    pub fn pearson_symbol(&self) -> String {
        self.0.pearson_symbol.clone()
    }

    // ------------------------------------------------------------------------
    // Primitive standardized layer cell
    // ------------------------------------------------------------------------
    /// Primitive standardized layer cell.
    #[getter]
    pub fn prim_std_cell(&self) -> PyStructure {
        layer_cell_to_cell(&self.0.prim_std_cell).into()
    }

    /// Linear part of the transformation from the input cell to the primitive
    /// standardized layer cell.
    #[getter]
    pub fn prim_std_linear(&self) -> [[f64; 3]; 3] {
        self.0.prim_std_linear_as_array()
    }

    /// Origin shift of the transformation from the input cell to the
    /// primitive standardized layer cell.
    #[getter]
    pub fn prim_std_origin_shift(&self) -> [f64; 3] {
        self.0.prim_std_origin_shift_as_array()
    }

    /// Mapping of sites in the input cell to those in the primitive
    /// standardized layer cell.
    #[getter]
    pub fn mapping_std_prim(&self) -> Vec<usize> {
        self.0.mapping_std_prim.clone()
    }

    // ------------------------------------------------------------------------
    // Final parameters
    // ------------------------------------------------------------------------
    /// Actual ``symprec`` used in the symmetry search.
    #[getter]
    pub fn symprec(&self) -> f64 {
        self.0.symprec
    }

    /// Actual ``angle_tolerance`` used in the symmetry search.
    #[getter]
    pub fn angle_tolerance(&self) -> Option<f64> {
        if let AngleTolerance::Radian(angle_tolerance) = self.0.angle_tolerance {
            Some(angle_tolerance)
        } else {
            None
        }
    }

    // ------------------------------------------------------------------------
    // Special methods
    // ------------------------------------------------------------------------
    fn __str__(&self) -> String {
        format!(
            "MoyoLayerDataset(number={}, hall_number={}, operations=<{} operations>, orbits={:?}, wyckoffs={:?}, site_symmetry_symbols={:?})",
            self.0.number,
            self.0.hall_number,
            self.0.operations.len(),
            self.0.orbits,
            self.0.wyckoffs,
            self.0.site_symmetry_symbols
        )
    }

    fn __repr__(&self) -> String {
        self.serialize_json()
    }

    // ------------------------------------------------------------------------
    // Serialization
    // ------------------------------------------------------------------------
    /// Serialize this object to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.0).expect("Serialization should not fail")
    }

    /// Deserialize an object from a JSON string.
    #[classmethod]
    pub fn deserialize_json(_cls: &Bound<'_, PyType>, s: &str) -> PyResult<Self> {
        serde_json::from_str(s).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Convert this object to a dictionary.
    pub fn as_dict(&self) -> PyResult<Py<PyAny>> {
        Python::attach(|py| {
            let obj = pythonize(py, &self.0).expect("Python object conversion should not fail");
            obj.into_py_any(py)
        })
    }

    /// Create an object from a dictionary.
    #[classmethod]
    pub fn from_dict(_cls: &Bound<'_, PyType>, obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Python::attach(|_| {
            depythonize::<Self>(obj).map_err(|e| {
                PyErr::new::<PyValueError, _>(format!("Deserialization failed: {}", e))
            })
        })
    }
}
