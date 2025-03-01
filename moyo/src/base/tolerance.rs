use log::debug;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

use super::error::MoyoError;

pub const EPS: f64 = 1e-8;

const INITIAL_SYMMETRY_SEARCH_STRIDE: f64 = 2.0;

#[derive(Debug, Copy, Clone, Deserialize, Serialize)]
/// Tolerance for angle in comparing basis vectors in symmetry search.
pub enum AngleTolerance {
    /// Tolerance in radian.
    Radian(f64),
    /// Default tolerance same as Spglib.
    Default,
}

pub trait Tolerances {
    fn increase_tolerances(&self, stride: f64) -> Self;
    fn reduce_tolerances(&self, stride: f64) -> Self;
}

#[derive(Debug, Clone)]
pub struct SymmetryTolerances {
    pub symprec: f64,
    pub angle_tolerance: AngleTolerance,
}

impl Tolerances for SymmetryTolerances {
    fn increase_tolerances(&self, stride: f64) -> Self {
        let symprec = self.symprec * stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle * stride)
        } else {
            AngleTolerance::Default
        };
        Self {
            symprec,
            angle_tolerance,
        }
    }

    fn reduce_tolerances(&self, stride: f64) -> Self {
        let symprec = self.symprec / stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle / stride)
        } else {
            AngleTolerance::Default
        };
        Self {
            symprec,
            angle_tolerance,
        }
    }
}

#[derive(Debug, Clone)]
pub struct MagneticSymmetryTolerances {
    pub symprec: f64,
    pub angle_tolerance: AngleTolerance,
    pub mag_symprec: f64,
}

impl Tolerances for MagneticSymmetryTolerances {
    fn increase_tolerances(&self, stride: f64) -> Self {
        let symprec = self.symprec * stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle * stride)
        } else {
            AngleTolerance::Default
        };
        let mag_symprec = self.mag_symprec * stride;
        Self {
            symprec,
            angle_tolerance,
            mag_symprec,
        }
    }

    fn reduce_tolerances(&self, stride: f64) -> Self {
        let symprec = self.symprec / stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle / stride)
        } else {
            AngleTolerance::Default
        };
        let mag_symprec = self.mag_symprec / stride;
        Self {
            symprec,
            angle_tolerance,
            mag_symprec,
        }
    }
}

pub struct ToleranceHandler<T: Tolerances> {
    pub tolerances: T,
    stride: f64,
    prev_error: Option<MoyoError>,
}

impl<T: Tolerances + Debug> ToleranceHandler<T> {
    pub fn new(tolerances: T) -> Self {
        Self {
            tolerances,
            stride: INITIAL_SYMMETRY_SEARCH_STRIDE,
            prev_error: None,
        }
    }

    pub fn update(&mut self, err: MoyoError) {
        // Update stride
        if self.prev_error.is_some() && self.prev_error != Some(err) {
            self.stride = self.stride.sqrt()
        }
        self.prev_error = Some(err);

        // Update tolerances
        self.tolerances = match err {
            MoyoError::TooSmallToleranceError => {
                let new_tolerances = self.tolerances.increase_tolerances(self.stride);
                debug!("Increase tolerances: {:?}", new_tolerances);
                new_tolerances
            }
            _ => {
                let new_tolerances = self.tolerances.reduce_tolerances(self.stride);
                debug!("Reduce tolerances: {:?}", new_tolerances);
                new_tolerances
            }
        }
    }
}
