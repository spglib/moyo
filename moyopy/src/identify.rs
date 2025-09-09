use pyo3::prelude::*;

use moyo::base::Lattice;
use moyo::data::ArithmeticNumber;
use moyo::identify::PointGroup;

use crate::base::PyMoyoError;
use crate::utils::{to_3x3_slice, to_matrix3};

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
