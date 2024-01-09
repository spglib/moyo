use nalgebra::base::{Matrix3, Vector3};

use super::cell::Cell;
use super::lattice::Lattice;
use super::operation::{Operations, Rotation, Translation};

pub type UnimodularLinear = Matrix3<i32>;
pub type Linear = Matrix3<f64>;
pub type OriginShift = Vector3<f64>;

/// Represent change of origin and basis for an affine space
#[derive(Debug, Clone)]
pub struct UnimodularTransformation {
    pub linear: UnimodularLinear,
    pub origin_shift: OriginShift,
    // Inverse of unimodular matrix is also unimodular
    linear_inv: UnimodularLinear,
}

impl UnimodularTransformation {
    pub fn new(linear: UnimodularLinear, origin_shift: OriginShift) -> Self {
        let size = linear.map(|e| e as f64).determinant().round().abs() as usize;
        if size != 1 {
            panic!("transformation matrix should be unimodular");
        }

        let linear_inv = linear
            .map(|e| e as f64)
            .try_inverse()
            .unwrap()
            .map(|e| e.round() as i32);

        Self {
            linear,
            origin_shift,
            linear_inv,
        }
    }

    pub fn from_linear(linear: UnimodularLinear) -> Self {
        Self::new(linear, OriginShift::zeros())
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }

    pub fn transform_lattice(&self, lattice: &Lattice) -> Lattice {
        Lattice::new(lattice.basis * self.linear_as_f64())
    }

    pub fn transform_operations(&self, operations: &Operations) -> Operations {
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in operations
            .rotations
            .iter()
            .zip(operations.translations.iter())
        {
            let new_rotation = self.linear_inv * rotation * self.linear;
            let new_translation = self.linear_inv.map(|e| e as f64)
                * (rotation.map(|e| e as f64) * self.origin_shift + translation
                    - self.origin_shift);
            new_rotations.push(new_rotation);
            new_translations.push(new_translation);
        }
        Operations::new(new_rotations, new_translations)
    }

    pub fn transform_cell(&self, cell: &Cell) -> Cell {
        let new_lattice = self.transform_lattice(&cell.lattice);
        // (P, p)^-1 x = P^-1 (x - p)
        let new_positions = cell
            .positions
            .iter()
            .map(|pos| self.linear_inv.map(|e| e as f64) * (pos - self.origin_shift))
            .collect();
        Cell::new(new_lattice, new_positions, cell.numbers.clone())
    }
}

/// Represent change of origin and basis for an affine space
#[derive(Debug, Clone)]
pub struct Transformation {
    pub linear: Linear,
    pub origin_shift: OriginShift,
    linear_inv: Linear,
}

impl Transformation {
    pub fn new(linear: Linear, origin_shift: OriginShift) -> Self {
        let linear_inv = linear.map(|e| e as f64).try_inverse().unwrap();

        Self {
            linear,
            origin_shift,
            linear_inv,
        }
    }

    pub fn from_linear(linear: Linear) -> Self {
        Self::new(linear, OriginShift::zeros())
    }

    pub fn transform_lattice(&self, lattice: &Lattice) -> Lattice {
        Lattice::new(lattice.basis * self.linear)
    }

    /// (P, p)^-1 (W, w) (P, p)
    /// This function may decrease the number of operations if the transformation is not compatible with an operation.
    pub fn transform_operations(&self, operations: &Operations) -> Operations {
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in operations
            .rotations
            .iter()
            .zip(operations.translations.iter())
        {
            if let Some((new_rotation, new_translation)) = transform_operation(
                rotation,
                translation,
                &self.linear,
                &self.linear_inv,
                &self.origin_shift,
            ) {
                new_rotations.push(new_rotation);
                new_translations.push(new_translation);
            }
        }
        Operations::new(new_rotations, new_translations)
    }

    // Apply `trans`, which may increase the number of atoms in the cell.
    // Mapping from sites of the new cell to those of the original cell is also returned.
    // pub fn expand_transform(&self, transformation: &Transformation) -> (Self, SiteMapping) {
    //     unimplemented!()
    // }
}

fn transform_operation(
    rotation: &Rotation,
    translation: &Translation,
    linear: &Linear,
    linear_inv: &Linear,
    origin_shift: &OriginShift,
) -> Option<(Rotation, Translation)> {
    let new_rotation = (linear_inv * rotation.map(|e| e as f64) * linear).map(|e| e.round() as i32);

    // Check if `new_rotation` is an integer matrix
    let recovered =
        (linear * new_rotation.map(|e| e as f64) * linear_inv).map(|e| e.round() as i32);
    if recovered != *rotation {
        return None;
    }

    let new_translation =
        linear_inv * (rotation.map(|e| e as f64) * origin_shift + translation - origin_shift);
    Some((new_rotation, new_translation))
}

#[cfg(test)]
mod tests {
    use std::vec;

    use nalgebra::matrix;

    use super::Transformation;
    use crate::base::operation::{Operations, Translation};

    #[test]
    fn test_incompatible_transformation() {
        let transformation = Transformation::from_linear(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 2.0;
        ]);
        // threefold rotation
        let operations = Operations::new(
            vec![matrix![
                0, 0, 1;
                1, 0, 0;
                0, 1, 0;
            ]],
            vec![Translation::zeros()],
        );
        assert_eq!(
            transformation
                .transform_operations(&operations)
                .num_operations(),
            0
        );
    }
}
