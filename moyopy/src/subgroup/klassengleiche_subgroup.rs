use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use moyo::base::Operation;
use moyo::subgroup::{
    KlassengleicheSubgroup, KlassengleicheSubgroupConjugacyClass, KlassengleicheSubgroupConjugate,
    enumerate_klassengleiche_subgroups as moyo_enumerate_klassengleiche_subgroups,
    enumerate_klassengleiche_subgroups_by_index as moyo_enumerate_klassengleiche_subgroups_by_index,
};
use moyo::utils::{to_3x3_slice, to_matrix3, to_vector3};

use crate::base::{PyMoyoError, PyOperations};

/// A Klassengleiche subgroup embedded in its parent space group.
#[derive(Debug, Clone)]
#[pyclass(name = "KlassengleicheSubgroup", frozen, skip_from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyKlassengleicheSubgroup(KlassengleicheSubgroup);

#[pymethods]
impl PyKlassengleicheSubgroup {
    /// Coset representatives in the subgroup primitive basis.
    #[getter]
    pub fn operations(&self) -> PyOperations {
        self.0.operations.clone().into()
    }

    /// The same symmetry operations expressed in the parent primitive basis.
    #[getter]
    pub fn parent_operations(&self) -> PyOperations {
        self.0.parent_operations.clone().into()
    }

    /// Mapping from ``operations`` to the input parent operations.
    #[getter]
    pub fn operation_indices(&self) -> Vec<usize> {
        self.0.operation_indices.clone()
    }

    /// Basis transformation whose columns generate the translation sublattice.
    #[getter]
    pub fn transformation(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.transformation)
    }

    /// Klassengleiche index, equal to the translation-sublattice index.
    #[getter]
    pub fn klassengleiche_index(&self) -> usize {
        self.0.klassengleiche_index
    }

    fn __repr__(&self) -> String {
        format!(
            "KlassengleicheSubgroup(klassengleiche_index={}, transformation={:?}, operations=<{} operations>)",
            self.0.klassengleiche_index,
            to_3x3_slice(&self.0.transformation),
            self.0.operations.len()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<KlassengleicheSubgroup> for PyKlassengleicheSubgroup {
    fn from(inner: KlassengleicheSubgroup) -> Self {
        Self(inner)
    }
}

/// A member of a parent-conjugacy class of embedded Klassengleiche subgroups.
#[derive(Debug, Clone)]
#[pyclass(name = "KlassengleicheSubgroupConjugate", frozen, skip_from_py_object)]
#[pyo3(module = "moyopy")]
pub struct PyKlassengleicheSubgroupConjugate(KlassengleicheSubgroupConjugate);

#[pymethods]
impl PyKlassengleicheSubgroupConjugate {
    /// The conjugate subgroup embedded in the input parent.
    #[getter]
    pub fn subgroup(&self) -> PyKlassengleicheSubgroup {
        self.0.subgroup.clone().into()
    }

    /// The parent-space-group operation whose inverse gives this conjugate.
    ///
    /// This ``Operations`` object contains exactly one operation.
    #[getter]
    pub fn conjugator(&self) -> PyOperations {
        vec![self.0.conjugator.clone()].into()
    }

    fn __repr__(&self) -> String {
        format!(
            "KlassengleicheSubgroupConjugate(subgroup={:?}, conjugator=<1 operation>)",
            PyKlassengleicheSubgroup::from(self.0.subgroup.clone())
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<KlassengleicheSubgroupConjugate> for PyKlassengleicheSubgroupConjugate {
    fn from(inner: KlassengleicheSubgroupConjugate) -> Self {
        Self(inner)
    }
}

/// A parent-conjugacy class of embedded Klassengleiche subgroups.
#[derive(Debug, Clone)]
#[pyclass(
    name = "KlassengleicheSubgroupConjugacyClass",
    frozen,
    skip_from_py_object
)]
#[pyo3(module = "moyopy")]
pub struct PyKlassengleicheSubgroupConjugacyClass(KlassengleicheSubgroupConjugacyClass);

#[pymethods]
impl PyKlassengleicheSubgroupConjugacyClass {
    /// Deterministically selected representative of this conjugacy class.
    #[getter]
    pub fn representative(&self) -> PyKlassengleicheSubgroup {
        self.0.representative.clone().into()
    }

    /// Every distinct conjugate embedding, including the representative first.
    #[getter]
    pub fn conjugates(&self) -> Vec<PyKlassengleicheSubgroupConjugate> {
        self.0
            .conjugates
            .iter()
            .cloned()
            .map(PyKlassengleicheSubgroupConjugate::from)
            .collect()
    }

    /// Coset representatives of the normalizer modulo the subgroup translation lattice.
    #[getter]
    pub fn normalizer_operations(&self) -> PyOperations {
        self.0.normalizer_operations.clone().into()
    }

    /// Permutations induced on the representative operations by the normalizer.
    ///
    /// This list is aligned with ``normalizer_operations``. In each permutation,
    /// operation position ``i`` is mapped to the value at position ``i``.
    #[getter]
    pub fn normalizer_permutations(&self) -> Vec<Vec<usize>> {
        self.0.normalizer_permutations.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "KlassengleicheSubgroupConjugacyClass(representative={:?}, conjugates=<{} subgroups>, normalizer=<{} operations>)",
            PyKlassengleicheSubgroup::from(self.0.representative.clone()),
            self.0.conjugates.len(),
            self.0.normalizer_operations.len()
        )
    }

    fn __str__(&self) -> String {
        self.__repr__()
    }
}

impl From<KlassengleicheSubgroupConjugacyClass> for PyKlassengleicheSubgroupConjugacyClass {
    fn from(inner: KlassengleicheSubgroupConjugacyClass) -> Self {
        Self(inner)
    }
}

fn operations_from_arrays(
    prim_rotations: &[[[i32; 3]; 3]],
    prim_translations: &[[f64; 3]],
) -> PyResult<Vec<Operation>> {
    if prim_rotations.len() != prim_translations.len() {
        return Err(PyValueError::new_err(
            "prim_rotations and prim_translations must have the same length",
        ));
    }

    Ok(prim_rotations
        .iter()
        .zip(prim_translations)
        .map(|(rotation, translation)| {
            Operation::new(to_matrix3(rotation), to_vector3(translation))
        })
        .collect())
}

/// Enumerate Klassengleiche subgroups with a prescribed translation sublattice.
///
/// ``transformation`` is expressed in the input primitive basis, and its columns
/// generate the subgroup translation lattice. It must have positive determinant
/// and be preserved by every parent rotation. The returned subgroups are
/// partitioned into conjugacy classes under the parent space group.
#[pyfunction]
#[pyo3(signature = (prim_rotations, prim_translations, transformation, *, epsilon=1e-4))]
pub fn enumerate_klassengleiche_subgroups(
    prim_rotations: Vec<[[i32; 3]; 3]>,
    prim_translations: Vec<[f64; 3]>,
    transformation: [[i32; 3]; 3],
    epsilon: f64,
) -> PyResult<Vec<PyKlassengleicheSubgroupConjugacyClass>> {
    let prim_operations = operations_from_arrays(&prim_rotations, &prim_translations)?;
    let transformation = to_matrix3(&transformation);
    let classes =
        moyo_enumerate_klassengleiche_subgroups(&prim_operations, &transformation, epsilon)
            .map_err(PyMoyoError::from)?
            .into_iter()
            .map(PyKlassengleicheSubgroupConjugacyClass::from)
            .collect();
    Ok(classes)
}

/// Enumerate Klassengleiche subgroups with a prescribed index.
///
/// Every index-``index`` translation sublattice is considered. Sublattices not
/// preserved by the parent point group are skipped, and the returned subgroups
/// retain their exact embedding in the input primitive basis.
#[pyfunction]
#[pyo3(signature = (prim_rotations, prim_translations, index, *, epsilon=1e-4))]
pub fn enumerate_klassengleiche_subgroups_by_index(
    prim_rotations: Vec<[[i32; 3]; 3]>,
    prim_translations: Vec<[f64; 3]>,
    index: usize,
    epsilon: f64,
) -> PyResult<Vec<PyKlassengleicheSubgroupConjugacyClass>> {
    let prim_operations = operations_from_arrays(&prim_rotations, &prim_translations)?;
    let classes =
        moyo_enumerate_klassengleiche_subgroups_by_index(&prim_operations, index, epsilon)
            .map_err(PyMoyoError::from)?
            .into_iter()
            .map(PyKlassengleicheSubgroupConjugacyClass::from)
            .collect();
    Ok(classes)
}
