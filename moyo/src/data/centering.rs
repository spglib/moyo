use nalgebra::{Matrix3, Vector3};
use strum_macros::EnumIter;

use crate::base::{Linear, Translation};

#[derive(Debug, Copy, Clone, PartialEq, EnumIter)]
pub enum Centering {
    P, // Primitive
    A, // A-face centered
    B, // B-face centered
    C, // C-face centered
    I, // Body centered
    R, // Rhombohedral (obverse setting)
    F, // Face centered
}

impl Centering {
    pub fn order(&self) -> usize {
        match self {
            Centering::P => 1,
            Centering::A => 2,
            Centering::B => 2,
            Centering::C => 2,
            Centering::I => 2,
            Centering::R => 3,
            Centering::F => 4,
        }
    }

    /// Transformation matrix from primitive to conventional cell.
    // Inverse matrices of https://github.com/spglib/spglib/blob/39a95560dd831c2d16f162126921ac1e519efa31/src/spacegroup.c#L373-L384
    pub fn linear(&self) -> Linear {
        match self {
            Centering::P => Linear::identity(),
            Centering::A => Linear::new(
                1, 0, 0, //
                0, 1, 1, //
                0, -1, 1, //
            ),
            Centering::B => Linear::new(
                1, 0, -1, //
                0, 1, 0, //
                1, 0, 1, //
            ),
            Centering::C => Linear::new(
                1, -1, 0, //
                1, 1, 0, //
                0, 0, 1, //
            ),
            Centering::R => Linear::new(
                1, 0, 1, //
                -1, 1, 1, //
                0, -1, 1, //
            ),
            Centering::I => Linear::new(
                0, 1, 1, //
                1, 0, 1, //
                1, 1, 0, //
            ),
            Centering::F => Linear::new(
                -1, 1, 1, //
                1, -1, 1, //
                1, 1, -1, //
            ),
        }
    }

    /// Transformation matrix from conventional to primitive cell.
    pub fn inverse(&self) -> Matrix3<f64> {
        self.linear().map(|e| e as f64).try_inverse().unwrap()
    }

    pub fn lattice_points(&self) -> Vec<Vector3<f64>> {
        match self {
            Centering::P => {
                vec![Translation::new(0.0, 0.0, 0.0)]
            }
            Centering::A => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                ]
            }
            Centering::B => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.0, 0.5),
                ]
            }
            Centering::C => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.0),
                ]
            }
            Centering::I => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.5, 0.5, 0.5),
                ]
            }
            Centering::R => {
                // obverse setting
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
                    Translation::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
                ]
            }
            // Centering::H => {
            //     vec![
            //         Translation::new(0.0, 0.0, 0.0),
            //         Translation::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
            //         Translation::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
            //     ]
            // }
            Centering::F => {
                vec![
                    Translation::new(0.0, 0.0, 0.0),
                    Translation::new(0.0, 0.5, 0.5),
                    Translation::new(0.5, 0.0, 0.5),
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
        for centering in Centering::iter() {
            assert_eq!(
                Transformation::from_linear(centering.linear()).size,
                centering.order()
            );
        }
    }
}
