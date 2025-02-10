use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use super::action::RotationMagneticMomentAction;
use super::cell::{AtomicSpecie, Cell, Position};
use super::lattice::Lattice;
use super::operation::{CartesianRotation, TimeReversal};

pub trait MagneticMoment: Sized + Clone {
    fn act_rotation(
        &self,
        cartesian_rotation: &CartesianRotation,
        action: RotationMagneticMomentAction,
    ) -> Self;

    fn act_time_reversal(&self, time_reversal: TimeReversal) -> Self;

    fn is_close(&self, other: &Self, mag_symprec: f64) -> bool;

    fn average(magnetic_moments: &[Self]) -> Self;

    fn act_magnetic_operation(
        &self,
        cartesian_rotation: &CartesianRotation,
        time_reversal: TimeReversal,
        action: RotationMagneticMomentAction,
    ) -> Self {
        let rotated = self.act_rotation(cartesian_rotation, action);
        rotated.act_time_reversal(time_reversal)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Collinear(pub f64);

impl MagneticMoment for Collinear {
    fn act_rotation(
        &self,
        cartesian_rotation: &CartesianRotation,
        action: RotationMagneticMomentAction,
    ) -> Self {
        match action {
            RotationMagneticMomentAction::Polar => Self(self.0),
            RotationMagneticMomentAction::Axial => {
                let det = cartesian_rotation.determinant().round();
                Self(det * self.0)
            }
        }
    }

    fn act_time_reversal(&self, time_reversal: TimeReversal) -> Self {
        if time_reversal {
            Self(-self.0)
        } else {
            Self(self.0)
        }
    }

    fn is_close(&self, other: &Self, mag_symprec: f64) -> bool {
        (self.0 - other.0).abs() < mag_symprec
    }

    fn average(magnetic_moments: &[Self]) -> Self {
        let sum = magnetic_moments.iter().map(|m| m.0).sum::<f64>();
        Collinear(sum / magnetic_moments.len() as f64)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NonCollinear(pub Vector3<f64>);

impl MagneticMoment for NonCollinear {
    fn act_rotation(
        &self,
        cartesian_rotation: &CartesianRotation,
        action: RotationMagneticMomentAction,
    ) -> Self {
        match action {
            RotationMagneticMomentAction::Polar => Self(cartesian_rotation * self.0),
            RotationMagneticMomentAction::Axial => {
                let det = cartesian_rotation.determinant().round();
                Self(det * cartesian_rotation * self.0)
            }
        }
    }

    fn act_time_reversal(&self, time_reversal: TimeReversal) -> Self {
        if time_reversal {
            Self(-self.0)
        } else {
            Self(self.0)
        }
    }

    fn is_close(&self, other: &Self, mag_symprec: f64) -> bool {
        // L2 norm
        (self.0 - other.0).norm() < mag_symprec
    }

    fn average(magnetic_moments: &[Self]) -> Self {
        let sum = magnetic_moments
            .iter()
            .map(|m| m.0)
            .fold(Vector3::zeros(), |acc, x| acc + x);
        NonCollinear(sum / magnetic_moments.len() as f64)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MagneticCell<M: MagneticMoment> {
    pub cell: Cell,
    pub magnetic_moments: Vec<M>,
}

impl<M: MagneticMoment> MagneticCell<M> {
    pub fn new(
        lattice: Lattice,
        positions: Vec<Position>,
        numbers: Vec<AtomicSpecie>,
        magnetic_moments: Vec<M>,
    ) -> Self {
        let cell = Cell::new(lattice, positions, numbers);
        Self::from_cell(cell, magnetic_moments)
    }

    pub fn from_cell(cell: Cell, magnetic_moments: Vec<M>) -> Self {
        if cell.positions.len() != magnetic_moments.len() {
            panic!("positions and magnetic_moments should be the same length");
        }

        Self {
            cell,
            magnetic_moments,
        }
    }

    pub fn num_atoms(&self) -> usize {
        self.cell.num_atoms()
    }
}
