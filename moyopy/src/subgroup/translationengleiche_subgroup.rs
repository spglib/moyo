use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::base::Operation;
use moyo::subgroup::{
    TranslationengleicheSubgroup, TranslationengleicheSubgroupConjugacyClass,
    TranslationengleicheSubgroupConjugate,
    enumerate_translationengleiche_subgroups as identify_enumerate_translationengleiche_subgroups,
};
use moyo::utils::{to_matrix3, to_vector3};

use crate::base::{PyMoyoError, PyOperations};

/// A Translationengleiche subgroup embedded in its parent space group.
#[derive(Debug, Clone)]
#[pyclass(name = "TranslationengleicheSubgroup", frozen, skip_from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyTranslationengleicheSubgroup(TranslationengleicheSubgroup);

#[pymethods]
impl PyTranslationengleicheSubgroup {
    /// Coset representatives of the subgroup, in the parent primitive basis.
    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

    /// Mapping from ``operations`` to the input parent operations.
    #[getter]
    pub fn operation_indices(&self) -> Vec<usize> {
        self.0.operation_indices.clone()
    }

    /// Translationengleiche index ``[G/T : H/T]``.
    #[getter]
    pub fn translationengleiche_index(&self) -> usize {
        self.0.translationengleiche_index
    }

    /// Shortest covering-chain length from the parent in the subgroup lattice.
    #[getter]
    pub fn depth(&self) -> usize {
        self.0.depth
    }

    fn __repr__(&self) -> String {
        format!(
            "TranslationengleicheSubgroup(translationengleiche_index={}, depth={}, operations=<{} operations>)",
            self.0.translationengleiche_index,
            self.0.depth,
            self.0.operations.len()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<TranslationengleicheSubgroup> for PyTranslationengleicheSubgroup {
    fn from(inner: TranslationengleicheSubgroup) -> Self {
        Self(inner)
    }
}

/// A member of a parent-conjugacy class of embedded Translationengleiche subgroups.
#[derive(Debug, Clone)]
#[pyclass(
    name = "TranslationengleicheSubgroupConjugate",
    frozen,
    skip_from_py_object
)]
#[pyo3(module = "moyopy")]
pub struct PyTranslationengleicheSubgroupConjugate(TranslationengleicheSubgroupConjugate);

#[pymethods]
impl PyTranslationengleicheSubgroupConjugate {
    /// The conjugate subgroup embedded in the input parent.
    #[getter]
    pub fn subgroup(&self) -> PyTranslationengleicheSubgroup {
        self.0.subgroup.clone().into()
    }

    /// Index of an input parent operation whose inverse conjugates the class
    /// representative to ``subgroup``.
    #[getter]
    pub fn conjugator_index(&self) -> usize {
        self.0.conjugator_index
    }

    fn __repr__(&self) -> String {
        format!(
            "TranslationengleicheSubgroupConjugate(conjugator_index={}, subgroup={:?})",
            self.0.conjugator_index,
            PyTranslationengleicheSubgroup::from(self.0.subgroup.clone())
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<TranslationengleicheSubgroupConjugate> for PyTranslationengleicheSubgroupConjugate {
    fn from(inner: TranslationengleicheSubgroupConjugate) -> Self {
        Self(inner)
    }
}

/// A parent-conjugacy class of embedded Translationengleiche subgroups.
#[derive(Debug, Clone)]
#[pyclass(
    name = "TranslationengleicheSubgroupConjugacyClass",
    frozen,
    skip_from_py_object
)]
#[pyo3(module = "moyopy")]
pub struct PyTranslationengleicheSubgroupConjugacyClass(TranslationengleicheSubgroupConjugacyClass);

#[pymethods]
impl PyTranslationengleicheSubgroupConjugacyClass {
    /// Deterministically selected representative of this conjugacy class.
    #[getter]
    pub fn representative(&self) -> PyTranslationengleicheSubgroup {
        self.0.representative.clone().into()
    }

    /// Every distinct conjugate embedding, including the representative first.
    #[getter]
    pub fn conjugates(&self) -> Vec<PyTranslationengleicheSubgroupConjugate> {
        self.0
            .conjugates
            .iter()
            .cloned()
            .map(PyTranslationengleicheSubgroupConjugate::from)
            .collect()
    }

    /// Indices of input parent operations that normalize the representative.
    #[getter]
    pub fn normalizer_indices(&self) -> Vec<usize> {
        self.0.normalizer_indices.clone()
    }

    /// Permutations induced on the representative operations by the normalizer.
    ///
    /// This list is aligned with ``normalizer_indices``. In each permutation,
    /// operation position ``i`` is mapped to the value at position ``i``.
    #[getter]
    pub fn normalizer_permutations(&self) -> Vec<Vec<usize>> {
        self.0.normalizer_permutations.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "TranslationengleicheSubgroupConjugacyClass(representative={:?}, conjugates=<{} subgroups>, normalizer=<{} operations>)",
            PyTranslationengleicheSubgroup::from(self.0.representative.clone()),
            self.0.conjugates.len(),
            self.0.normalizer_indices.len()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<TranslationengleicheSubgroupConjugacyClass>
    for PyTranslationengleicheSubgroupConjugacyClass
{
    fn from(inner: TranslationengleicheSubgroupConjugacyClass) -> Self {
        Self(inner)
    }
}

/// Enumerate Translationengleiche subgroups of a primitive space group.
///
/// The result contains every Translationengleiche subgroup, including the
/// parent and pure-translation subgroup, partitioned into conjugacy classes
/// under the parent. Returned operations stay in the input primitive basis,
/// and operation and conjugator indices refer directly to the input arrays.
///
/// Parameters
/// ----------
/// prim_rotations : list\[list\[list\[int\]\]\]
///     One rotation matrix for every point-group element of the primitive cell.
/// prim_translations : list\[list\[float\]\]
///     Corresponding translation vectors in the primitive basis.
/// epsilon : float
///     Numerical tolerance for checking closure modulo lattice translations.
#[pyfunction]
#[pyo3(signature = (prim_rotations, prim_translations, *, epsilon=1e-4))]
pub fn enumerate_translationengleiche_subgroups(
    prim_rotations: Vec<[[i32; 3]; 3]>,
    prim_translations: Vec<[f64; 3]>,
    epsilon: f64,
) -> PyResult<Vec<PyTranslationengleicheSubgroupConjugacyClass>> {
    if prim_rotations.len() != prim_translations.len() {
        return Err(PyValueError::new_err(
            "prim_rotations and prim_translations must have the same length",
        ));
    }

    let prim_operations = prim_rotations
        .iter()
        .zip(&prim_translations)
        .map(|(rotation, translation)| {
            Operation::new(to_matrix3(rotation), to_vector3(translation))
        })
        .collect();
    let classes = identify_enumerate_translationengleiche_subgroups(&prim_operations, epsilon)
        .map_err(PyMoyoError::from)?
        .into_iter()
        .map(PyTranslationengleicheSubgroupConjugacyClass::from)
        .collect();
    Ok(classes)
}
