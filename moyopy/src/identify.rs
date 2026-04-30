use itertools::izip;
use pyo3::exceptions::PyValueError;
use pyo3::{IntoPyObjectExt, prelude::*};
use pythonize::pythonize;
use serde::Serialize;
use serde_json;

use moyo::base::{Lattice, MagneticOperation, Operation};
use moyo::data::{ArithmeticNumber, HallNumber, Number, Setting, UNINumber};
use moyo::identify::{
    MagneticSpaceGroup, PointGroup, SpaceGroup, integral_normalizer as identify_integral_normalizer,
};
use moyo::utils::{to_3x3_slice, to_matrix3, to_vector3};

use crate::base::{PyMoyoError, PyUnimodularTransformation};
use crate::data::PySetting;

fn has_same_rotation(lhs: &Operation, rhs: &Operation) -> bool {
    lhs.rotation == rhs.rotation
}

fn push_operation(operations: &mut Vec<Operation>, operation: Operation) -> bool {
    if operations
        .iter()
        .any(|candidate| has_same_rotation(candidate, &operation))
    {
        false
    } else {
        operations.push(operation);
        true
    }
}

fn generated_closure(generators: &[Operation], prim_operations: &[Operation]) -> Vec<Operation> {
    let mut closure = vec![];
    let mut cursor = 0;
    push_operation(&mut closure, Operation::identity());

    while cursor < closure.len() {
        let lhs = closure[cursor].clone();
        cursor += 1;

        for rhs in generators {
            let new_operation = lhs.clone() * rhs.clone();
            if prim_operations
                .iter()
                .any(|operation| has_same_rotation(operation, &new_operation))
            {
                push_operation(&mut closure, new_operation);
            }
        }
    }

    closure
}

fn derive_small_generators(prim_operations: &[Operation]) -> PyResult<Vec<Operation>> {
    let mut generators = vec![];
    let mut closure = generated_closure(&generators, prim_operations);

    for operation in prim_operations {
        if closure
            .iter()
            .any(|generated| has_same_rotation(generated, operation))
        {
            continue;
        }
        generators.push(operation.clone());
        closure = generated_closure(&generators, prim_operations);
    }

    if closure.len() != prim_operations.len() {
        return Err(PyValueError::new_err(
            "failed to derive generators whose closure matches prim_operations",
        ));
    }

    Ok(generators)
}

/// Compute the integral normalizer of a space group.
///
/// Returns the unimodular transformations that conjugate the given space group into itself.
///
/// Parameters
/// ----------
/// prim_rotations : list[list[list[int]]]
///     Rotation matrices of the symmetry operations of the primitive cell.
/// prim_translations : list[list[float]]
///     Translation vectors of the symmetry operations of the primitive cell.
/// prim_generators : list[int] | None
///     Optional indices of operations to use as generators. If ``None``, a small generating
///     set is derived automatically.
/// epsilon : float
///     Numerical tolerance for matching translations.
#[pyfunction]
#[pyo3(signature = (prim_rotations, prim_translations, *, prim_generators=None, epsilon=1e-4))]
pub fn integral_normalizer(
    prim_rotations: Vec<[[i32; 3]; 3]>,
    prim_translations: Vec<[f64; 3]>,
    prim_generators: Option<Vec<usize>>,
    epsilon: f64,
) -> PyResult<Vec<PyUnimodularTransformation>> {
    if prim_rotations.len() != prim_translations.len() {
        return Err(PyValueError::new_err(
            "prim_rotations and prim_translations must have the same length",
        ));
    }

    let prim_operations = prim_rotations
        .iter()
        .zip(prim_translations.iter())
        .map(|(rotation, translation)| {
            Operation::new(to_matrix3(rotation), to_vector3(translation))
        })
        .collect::<Vec<_>>();

    let prim_generators = if let Some(indices) = prim_generators {
        let mut generators = Vec::with_capacity(indices.len());
        for index in indices {
            let operation = prim_operations.get(index).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "prim_generators index {index} is out of range for {} operations",
                    prim_operations.len()
                ))
            })?;
            generators.push(operation.clone());
        }
        generators
    } else {
        derive_small_generators(&prim_operations)?
    };

    Ok(
        identify_integral_normalizer(&prim_operations, &prim_generators, epsilon)
            .into_iter()
            .map(PyUnimodularTransformation::from)
            .collect(),
    )
}

/// Point group identified from a list of rotation matrices.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "PointGroup", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyPointGroup(PointGroup);

#[pymethods]
impl PyPointGroup {
    /// Identify the point group of the given primitive rotations.
    ///
    /// Parameters
    /// ----------
    /// prim_rotations : list[list[list[int]]]
    ///     Rotation matrices in the primitive cell.
    /// basis : list[list[float]] | None
    ///     Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
    ///     is assumed.
    #[new]
    #[pyo3(signature = (prim_rotations, *, basis=None))]
    pub fn new(
        prim_rotations: Vec<[[i32; 3]; 3]>,
        basis: Option<[[f64; 3]; 3]>,
    ) -> Result<Self, PyMoyoError> {
        let prim_rotations = prim_rotations.iter().map(to_matrix3).collect::<Vec<_>>();
        let point_group = if let Some(basis) = basis {
            let lattice = Lattice::from_basis(basis);
            PointGroup::from_lattice(&lattice, &prim_rotations)?
        } else {
            PointGroup::new(&prim_rotations)?
        };

        Ok(Self(point_group))
    }

    /// Number for the arithmetic crystal class (1 - 73).
    #[getter]
    pub fn arithmetic_number(&self) -> ArithmeticNumber {
        self.0.arithmetic_number
    }

    /// Transformation matrix from the input primitive basis to the standardized basis.
    #[getter]
    pub fn prim_trans_mat(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.prim_trans_mat)
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

impl From<PyPointGroup> for PointGroup {
    fn from(point_group: PyPointGroup) -> Self {
        point_group.0
    }
}

impl From<PointGroup> for PyPointGroup {
    fn from(point_group: PointGroup) -> Self {
        PyPointGroup(point_group)
    }
}

/// Space group identified from a list of primitive symmetry operations.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "SpaceGroup", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PySpaceGroup(SpaceGroup);

#[pymethods]
impl PySpaceGroup {
    /// Identify the space group from primitive rotations and translations.
    ///
    /// Parameters
    /// ----------
    /// prim_rotations : list[list[list[int]]]
    ///     Rotation matrices of the symmetry operations of the primitive cell.
    /// prim_translations : list[list[float]]
    ///     Translation vectors of the symmetry operations of the primitive cell.
    /// basis : list[list[float]] | None
    ///     Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
    ///     is assumed.
    /// setting : Setting | None
    ///     Preference for the standardized setting of the detected space-group type.
    /// epsilon : float
    ///     Numerical tolerance for matching translations.
    #[new]
    #[pyo3(signature = (prim_rotations, prim_translations, *, basis=None, setting=None, epsilon=1e-4))]
    pub fn new(
        prim_rotations: Vec<[[i32; 3]; 3]>,
        prim_translations: Vec<[f64; 3]>,
        basis: Option<[[f64; 3]; 3]>,
        setting: Option<PySetting>,
        epsilon: f64,
    ) -> Result<Self, PyMoyoError> {
        let prim_operations = prim_rotations
            .iter()
            .zip(prim_translations.iter())
            .map(|(r, t)| Operation::new(to_matrix3(r), to_vector3(t)))
            .collect::<Vec<_>>();
        let setting = if let Some(setting) = setting {
            setting.into()
        } else {
            Setting::default()
        };

        let space_group = if let Some(basis) = basis {
            let lattice = Lattice::from_basis(basis);
            SpaceGroup::from_lattice(&lattice, &prim_operations, setting, epsilon)?
        } else {
            SpaceGroup::new(&prim_operations, setting, epsilon)?
        };

        Ok(Self(space_group))
    }

    /// ITA number for the identified space group (1 - 230).
    #[getter]
    pub fn number(&self) -> Number {
        self.0.number
    }

    /// Hall symbol number (1 - 530) for the chosen setting.
    #[getter]
    pub fn hall_number(&self) -> HallNumber {
        self.0.hall_number
    }

    /// Linear part of the transformation from the input primitive basis to the standardized
    /// basis.
    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        self.0.transformation.linear_as_array()
    }

    /// Origin shift of the transformation from the input primitive basis to the standardized
    /// basis.
    #[getter]
    pub fn origin_shift(&self) -> [f64; 3] {
        self.0.transformation.origin_shift_as_array()
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

impl From<PySpaceGroup> for SpaceGroup {
    fn from(space_group: PySpaceGroup) -> Self {
        space_group.0
    }
}

impl From<SpaceGroup> for PySpaceGroup {
    fn from(space_group: SpaceGroup) -> Self {
        PySpaceGroup(space_group)
    }
}

/// Magnetic space group identified from a list of primitive magnetic operations.
#[derive(Debug, Clone, Serialize)]
#[pyclass(name = "MagneticSpaceGroup", frozen, from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyMagneticSpaceGroup(MagneticSpaceGroup);

#[pymethods]
impl PyMagneticSpaceGroup {
    /// Identify the magnetic space group from primitive magnetic operations.
    ///
    /// Parameters
    /// ----------
    /// prim_rotations : list[list[list[int]]]
    ///     Rotation matrices of the magnetic operations of the primitive cell.
    /// prim_translations : list[list[float]]
    ///     Translation vectors of the magnetic operations of the primitive cell.
    /// prim_time_reversals : list[bool]
    ///     Time-reversal flag for each magnetic operation of the primitive cell.
    /// basis : list[list[float]] | None
    ///     Row-wise basis vectors of the primitive lattice. If ``None``, an identity basis
    ///     is assumed.
    /// epsilon : float
    ///     Numerical tolerance for matching translations.
    #[new]
    #[pyo3(signature = (prim_rotations, prim_translations, prim_time_reversals, *, basis=None, epsilon=1e-4))]
    pub fn new(
        prim_rotations: Vec<[[i32; 3]; 3]>,
        prim_translations: Vec<[f64; 3]>,
        prim_time_reversals: Vec<bool>,
        basis: Option<[[f64; 3]; 3]>,
        epsilon: f64,
    ) -> Result<Self, PyMoyoError> {
        let prim_mag_operations = izip!(
            prim_rotations.iter(),
            prim_translations.iter(),
            prim_time_reversals.iter(),
        )
        .map(|(rot, trans, tr)| MagneticOperation::new(to_matrix3(rot), to_vector3(trans), *tr))
        .collect::<Vec<_>>();

        let magnetic_space_group = if let Some(basis) = basis {
            let lattice = Lattice::from_basis(basis);
            MagneticSpaceGroup::from_lattice(&lattice, &prim_mag_operations, epsilon)?
        } else {
            MagneticSpaceGroup::new(&prim_mag_operations, epsilon)?
        };

        Ok(Self(magnetic_space_group))
    }

    /// Serial number of UNI (and BNS) symbols.
    #[getter]
    pub fn uni_number(&self) -> UNINumber {
        self.0.uni_number
    }

    /// Linear part of the transformation from the input primitive basis to the standardized
    /// basis.
    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        self.0.transformation.linear_as_array()
    }

    /// Origin shift of the transformation from the input primitive basis to the standardized
    /// basis.
    #[getter]
    pub fn origin_shift(&self) -> [f64; 3] {
        self.0.transformation.origin_shift_as_array()
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

impl From<PyMagneticSpaceGroup> for MagneticSpaceGroup {
    fn from(magnetic_space_group: PyMagneticSpaceGroup) -> Self {
        magnetic_space_group.0
    }
}

#[cfg(test)]
mod tests {
    use super::{derive_small_generators, generated_closure};
    use crate::data::operations_from_number;
    use moyo::base::Operations;

    fn primitive_operations_from_number(number: i32) -> Operations {
        operations_from_number(number, None, true).unwrap().into()
    }

    #[test]
    fn test_derive_small_generators_reaches_full_group() {
        let prim_operations = primitive_operations_from_number(221);
        let generators = derive_small_generators(&prim_operations).unwrap();
        let closure = generated_closure(&generators, &prim_operations);

        assert!(generators.len() < prim_operations.len());
        assert_eq!(closure.len(), prim_operations.len());
    }
}

impl From<MagneticSpaceGroup> for PyMagneticSpaceGroup {
    fn from(magnetic_space_group: MagneticSpaceGroup) -> Self {
        PyMagneticSpaceGroup(magnetic_space_group)
    }
}
