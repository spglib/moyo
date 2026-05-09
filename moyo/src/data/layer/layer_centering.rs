use nalgebra::{Matrix3, Vector3};
use serde::Serialize;
use strum_macros::EnumIter;

use crate::base::{Linear, Translation};

/// Centering type for layer groups (only oblique-/rectangular-primitive and
/// rectangular-centered exist in 2D).
#[derive(Debug, Copy, Clone, PartialEq, EnumIter, Serialize)]
pub enum LayerCentering {
    P, // Primitive
    C, // C-face centered (in-plane only)
}

impl LayerCentering {
    pub fn order(&self) -> usize {
        match self {
            LayerCentering::P => 1,
            LayerCentering::C => 2,
        }
    }

    /// Transformation matrix from primitive to conventional cell.
    /// The c (aperiodic) axis is left untouched.
    pub fn linear(&self) -> Linear {
        match self {
            LayerCentering::P => Linear::identity(),
            LayerCentering::C => Linear::new(
                1, -1, 0, //
                1, 1, 0, //
                0, 0, 1, //
            ),
        }
    }

    /// Transformation matrix from conventional to primitive cell.
    pub fn inverse(&self) -> Matrix3<f64> {
        self.linear().map(|e| e as f64).try_inverse().unwrap()
    }

    /// Lattice points of the conventional cell expressed in the conventional basis.
    /// Centerings are 2D-only: the c-component is always zero.
    pub fn lattice_points(&self) -> Vec<Vector3<f64>> {
        match self {
            LayerCentering::P => {
                vec![Translation::new(0.0, 0.0, 0.0)]
            }
            LayerCentering::C => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use strum::IntoEnumIterator;

    use super::*;
    use crate::base::Transformation;

    #[test]
    fn test_conventional_transformation_matrix() {
        for centering in LayerCentering::iter() {
            assert_eq!(
                Transformation::from_linear(centering.linear()).size,
                centering.order()
            );
        }
    }

    #[test]
    fn test_lattice_points_in_plane() {
        for centering in LayerCentering::iter() {
            for p in centering.lattice_points() {
                assert_eq!(p[2], 0.0);
            }
            assert_eq!(centering.lattice_points().len(), centering.order());
        }
    }
}
