use itertools::iproduct;
use nalgebra::base::{Matrix3, Vector3};

use super::cell::Cell;
use super::lattice::Lattice;
use super::operation::{Operations, Rotation, Translation};
use crate::math::SNF;

pub type UnimodularLinear = Matrix3<i32>;
pub type Linear = Matrix3<i32>;
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
        let det = linear.map(|e| e as f64).determinant().round() as i32;
        if det != 1 {
            panic!("Determinant of transformation matrix should be one.");
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

    pub fn from_origin_shift(origin_shift: OriginShift) -> Self {
        Self::new(UnimodularLinear::identity(), origin_shift)
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }

    pub fn transform_lattice(&self, lattice: &Lattice) -> Lattice {
        Lattice::new((lattice.basis * self.linear_as_f64()).transpose())
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
    pub size: usize,
    pub linear_inv: Matrix3<f64>,
}

impl Transformation {
    pub fn new(linear: Linear, origin_shift: OriginShift) -> Self {
        let linear_inv = linear.map(|e| e as f64).try_inverse().unwrap();

        let det = linear.map(|e| e as f64).determinant().round() as i32;
        if det <= 0 {
            panic!("Determinant of transformation matrix should be positive.");
        }

        Self {
            linear,
            origin_shift,
            size: det as usize,
            linear_inv,
        }
    }

    pub fn from_linear(linear: Linear) -> Self {
        Self::new(linear, OriginShift::zeros())
    }

    pub fn from_origin_shift(origin_shift: OriginShift) -> Self {
        Self::new(Linear::identity(), origin_shift)
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }

    pub fn transform_lattice(&self, lattice: &Lattice) -> Lattice {
        self.transform_lattice_with_linear(lattice, &self.linear_as_f64())
    }

    pub fn inverse_transform_lattice(&self, lattice: &Lattice) -> Lattice {
        self.transform_lattice_with_linear(lattice, &self.linear_inv)
    }

    fn transform_lattice_with_linear(&self, lattice: &Lattice, linear: &Matrix3<f64>) -> Lattice {
        Lattice::new((lattice.basis * linear).transpose())
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
            if let Some((new_rotation, new_translation)) = transform_operation_as_f64(
                rotation,
                translation,
                &self.linear.map(|e| e as f64),
                &self.linear_inv,
                &self.origin_shift,
            ) {
                new_rotations.push(new_rotation);
                new_translations.push(new_translation);
            }
        }
        Operations::new(new_rotations, new_translations)
    }

    /// (P, p) (W, w) (P, p)^-1
    /// This function may decrease the number of operations if the transformation is not compatible with an operation.
    pub fn inverse_transform_operations(&self, operations: &Operations) -> Operations {
        let mut new_rotations = vec![];
        let mut new_translations = vec![];
        for (rotation, translation) in operations
            .rotations
            .iter()
            .zip(operations.translations.iter())
        {
            if let Some((new_rotation, new_translation)) = transform_operation_as_f64(
                rotation,
                translation,
                &self.linear_inv,
                &self.linear.map(|e| e as f64),
                &(-self.linear_as_f64() * self.origin_shift),
            ) {
                new_rotations.push(new_rotation);
                new_translations.push(new_translation);
            }
        }
        Operations::new(new_rotations, new_translations)
    }

    // The transformation may increase the number of atoms in the cell.
    // Return the transformed cell and mapping from sites in the transformed cell to sites in the original cell.
    pub fn transform_cell(&self, cell: &Cell) -> (Cell, Vec<usize>) {
        let new_lattice = self.transform_lattice(&cell.lattice);

        // Let self.linear be M.
        // Consider the Smith normal form of M: D = L * M * R. S = diag([D_0, D_1, D_2])
        // For `x` in unit cell and `n` in Z^3,
        //      A * (x + n) = (A * M) * M^-1 * (x + n)
        // Two lattice points with n - n' in M * Z^3 are equivalent:
        //      n = n' (mod M * Z^3) <=> M^-1 * (n - n') = 0 (mod Z^3)
        //                           <=> R * D^-1 * L * (n - n') = 0 (mod Z^3)
        //                           <=> L * n = L * n' (mod D * Z^3)
        // Thus, distinct lattice points in the sublattice are generated by
        //      { n = L^-1 * f | f in Z_{D_0} otimes Z_{D_1} otimes Z_{D_2} }.
        let snf = SNF::new(&self.linear);
        let linv = snf
            .l
            .map(|e| e as f64)
            .try_inverse()
            .unwrap()
            .map(|e| e.round() as i32);
        let lattice_points = iproduct!((0..snf.d[(0, 0)]), (0..snf.d[(1, 1)]), (0..snf.d[(2, 2)]))
            .map(|(f0, f1, f2)| (linv * Vector3::new(f0, f1, f2)).map(|e| e as f64))
            .collect::<Vec<_>>();

        let mut new_positions = vec![];
        let mut new_numbers = vec![];
        let mut site_mapping = vec![];
        for (i, (pos, number)) in cell.positions.iter().zip(cell.numbers.iter()).enumerate() {
            for lattice_point in lattice_points.iter() {
                // Fractional coordinates in the new sublattice
                let new_position = (self.linear_inv * (pos + lattice_point)).map(|e| e % 1.);
                new_positions.push(new_position);
                new_numbers.push(*number);
                site_mapping.push(i);
            }
        }

        (
            Cell::new(new_lattice, new_positions, new_numbers),
            site_mapping,
        )
    }
}

/// Transform operation (rotation, translation) by transformation (linear, origin_shift).
fn transform_operation_as_f64(
    rotation: &Rotation,
    translation: &Translation,
    linear: &Matrix3<f64>,
    linear_inv: &Matrix3<f64>,
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
            1, 0, 0;
            0, 1, 0;
            0, 0, 2;
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
