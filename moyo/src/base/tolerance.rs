use log::debug;

use super::error::MoyoError;

pub const EPS: f64 = 1e-8;

const INITIAL_SYMMETRY_SEARCH_STRIDE: f64 = 2.0;

#[derive(Debug, Copy, Clone)]
/// Tolerance for angle in comparing basis vectors in symmetry search.
pub enum AngleTolerance {
    /// Tolerance in radian.
    Radian(f64),
    /// Default tolerance same as Spglib.
    Default,
}

pub struct ToleranceHandler {
    pub symprec: f64,
    pub angle_tolerance: AngleTolerance,
    stride: f64,
    prev_error: Option<MoyoError>,
}

impl ToleranceHandler {
    pub fn new(symprec: f64, angle_tolerance: AngleTolerance) -> Self {
        Self {
            symprec,
            angle_tolerance,
            stride: INITIAL_SYMMETRY_SEARCH_STRIDE,
            prev_error: None,
        }
    }

    pub fn update(&mut self, err: MoyoError) {
        // Update stride
        if !self.prev_error.is_none() && self.prev_error != Some(err) {
            self.stride = self.stride.sqrt()
        }
        self.prev_error = Some(err);

        // Update tolerances
        (self.symprec, self.angle_tolerance) = match err {
            MoyoError::TooSmallToleranceError => self.increase_tolerance(),
            MoyoError::TooLargeToleranceError => self.reduce_tolerance(),
            _ => self.reduce_tolerance(),
        }
    }

    fn increase_tolerance(&self) -> (f64, AngleTolerance) {
        let symprec = self.symprec * self.stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle * self.stride)
        } else {
            AngleTolerance::Default
        };
        debug!(
            "Increase tolerance to symprec={}, angle_tolerance={:?}",
            symprec, angle_tolerance
        );
        (symprec, angle_tolerance)
    }

    fn reduce_tolerance(&self) -> (f64, AngleTolerance) {
        let symprec = self.symprec / self.stride;
        let angle_tolerance = if let AngleTolerance::Radian(angle) = self.angle_tolerance {
            AngleTolerance::Radian(angle / self.stride)
        } else {
            AngleTolerance::Default
        };
        debug!(
            "Reduce tolerance to symprec={}, angle_tolerance={:?}",
            symprec, angle_tolerance
        );
        (symprec, angle_tolerance)
    }
}
