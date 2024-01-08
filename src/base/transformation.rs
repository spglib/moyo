use nalgebra::base::{Matrix3, Vector3};

pub type UnimodularLinear = Matrix3<i32>;
pub type Linear = Matrix3<f64>;
pub type OriginShift = Vector3<f64>;

/// Represent change of origin and basis for an affine space
#[derive(Debug, Clone)]
pub struct UnimodularTransformation {
    pub linear: UnimodularLinear,
    pub origin_shift: OriginShift,
}

impl UnimodularTransformation {
    pub fn new(linear: UnimodularLinear, origin_shift: OriginShift) -> Self {
        let size = linear.map(|e| e as f64).determinant().round().abs() as usize;
        if size != 1 {
            panic!("transformation matrix should be unimodular");
        }

        Self {
            linear,
            origin_shift,
        }
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }
}

/// Represent change of origin and basis for an affine space
#[derive(Debug, Clone)]
pub struct Transformation {
    pub linear: Linear,
    pub origin_shift: OriginShift,
}

impl Transformation {
    pub fn new(linear: Linear, origin_shift: OriginShift) -> Self {
        Self {
            linear,
            origin_shift,
        }
    }
}
