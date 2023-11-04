use nalgebra::base::{Matrix3, Vector3};
use std::ops::Mul;

pub type TransformationMatrix = Matrix3<i32>;
pub type OriginShift = Vector3<f64>;

/// Represent change of origin and basis for an affine space
#[derive(Debug)]
pub struct Transformation {
    pub trans_mat: TransformationMatrix,
    pub origin_shift: OriginShift,
    pub size: usize,
}

impl Transformation {
    pub fn new(trans_mat: TransformationMatrix, origin_shift: OriginShift) -> Self {
        let size = trans_mat.map(|e| e as f64).determinant().round().abs() as usize;
        if size == 0 {
            panic!("transformation matrix should have nonzero determinant");
        }

        Self {
            trans_mat,
            origin_shift,
            size,
        }
    }

    pub fn identity() -> Self {
        Self::new(TransformationMatrix::identity(), OriginShift::zeros())
    }

    pub fn trans_mat_as_f64(&self) -> Matrix3<f64> {
        self.trans_mat.map(|e| e as f64)
    }
}

impl Mul for Transformation {
    type Output = Self;

    // Required method
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.trans_mat * rhs.trans_mat,
            self.trans_mat_as_f64() * rhs.origin_shift + self.origin_shift,
        )
    }
}
