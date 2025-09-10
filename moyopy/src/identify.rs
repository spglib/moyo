use pyo3::prelude::*;

use moyo::base::{Lattice, Operation};
use moyo::data::{ArithmeticNumber, HallNumber, Number, Setting};
use moyo::identify::{PointGroup, SpaceGroup};

use crate::base::PyMoyoError;
use crate::data::PySetting;
use crate::utils::{to_3_slice, to_3x3_slice, to_matrix3, to_vector3};

#[derive(Debug, Clone)]
#[pyclass(name = "PointGroup", frozen)]
#[pyo3(module = "moyopy")]
pub struct PyPointGroup(PointGroup);

#[pymethods]
impl PyPointGroup {
    #[new]
    #[pyo3(signature = (prim_rotations, *, basis=None))]
    pub fn new(
        prim_rotations: Vec<[[i32; 3]; 3]>,
        basis: Option<[[f64; 3]; 3]>,
    ) -> Result<Self, PyMoyoError> {
        let prim_rotations = prim_rotations
            .iter()
            .map(|x| to_matrix3(x))
            .collect::<Vec<_>>();
        let point_group = if let Some(basis) = basis {
            let lattice = Lattice::from_basis(basis);
            PointGroup::from_lattice(&lattice, &prim_rotations)?
        } else {
            PointGroup::new(&prim_rotations)?
        };

        Ok(Self(point_group))
    }

    #[getter]
    pub fn arithmetic_number(&self) -> ArithmeticNumber {
        self.0.arithmetic_number
    }

    #[getter]
    pub fn prim_trans_mat(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.prim_trans_mat)
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

#[derive(Debug, Clone)]
#[pyclass(name = "SpaceGroup", frozen)]
#[pyo3(module = "moyopy")]
pub struct PySpaceGroup(SpaceGroup);

#[pymethods]
impl PySpaceGroup {
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
            Setting::Spglib
        };

        let space_group = if let Some(basis) = basis {
            let lattice = Lattice::from_basis(basis);
            SpaceGroup::from_lattice(&lattice, &prim_operations, setting, epsilon)?
        } else {
            SpaceGroup::new(&prim_operations, setting, epsilon)?
        };

        Ok(Self(space_group))
    }

    #[getter]
    pub fn number(&self) -> Number {
        self.0.number
    }

    #[getter]
    pub fn hall_number(&self) -> HallNumber {
        self.0.hall_number
    }

    #[getter]
    pub fn linear(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.0.transformation.linear)
    }

    #[getter]
    pub fn origin_shift(&self) -> [f64; 3] {
        to_3_slice(&self.0.transformation.origin_shift)
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
