use std::ops::Mul;

use itertools::iproduct;
use nalgebra::base::{Matrix3, Vector3};

use super::cell::Cell;
use super::lattice::Lattice;
use super::magnetic_cell::{MagneticCell, MagneticMoment};
use super::operation::{MagneticOperation, Operation};
use crate::math::SNF;

pub type UnimodularLinear = Matrix3<i32>;
pub type Linear = Matrix3<i32>;
/// Origin shift in a crystallographic basis
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

    #[allow(dead_code)]
    pub fn from_origin_shift(origin_shift: OriginShift) -> Self {
        Self::new(UnimodularLinear::identity(), origin_shift)
    }

    pub fn inverse(&self) -> Self {
        // (P, p)^-1 = (P^-1, -P^-1 p)
        Self::new(
            self.linear_inv,
            -self.linear_inv.map(|e| e as f64) * self.origin_shift,
        )
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }

    pub fn transform_lattice(&self, lattice: &Lattice) -> Lattice {
        Lattice::new((lattice.basis * self.linear_as_f64()).transpose())
    }

    pub fn transform_operation(&self, operation: &Operation) -> Operation {
        let new_rotation = self.linear_inv * operation.rotation * self.linear;
        let new_translation = self.linear_inv.map(|e| e as f64)
            * (operation.rotation.map(|e| e as f64) * self.origin_shift + operation.translation
                - self.origin_shift);
        Operation::new(new_rotation, new_translation)
    }

    pub fn transform_operations(&self, operations: &[Operation]) -> Vec<Operation> {
        operations
            .iter()
            .map(|ops| self.transform_operation(ops))
            .collect()
    }

    pub fn transform_magnetic_operation(
        &self,
        magnetic_operation: &MagneticOperation,
    ) -> MagneticOperation {
        let new_operation = self.transform_operation(&magnetic_operation.operation);
        MagneticOperation::from_operation(new_operation, magnetic_operation.time_reversal)
    }

    pub fn transform_magnetic_operations(
        &self,
        magnetic_operations: &[MagneticOperation],
    ) -> Vec<MagneticOperation> {
        magnetic_operations
            .iter()
            .map(|mops| self.transform_magnetic_operation(mops))
            .collect()
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

    pub fn transform_magnetic_moments<M: MagneticMoment>(
        &self,
        magnetic_moments: &Vec<M>,
    ) -> Vec<M> {
        // Magnetic moments are not transformed
        magnetic_moments.clone()
    }

    pub fn transform_magnetic_cell<M: MagneticMoment>(
        &self,
        magnetic_cell: &MagneticCell<M>,
    ) -> MagneticCell<M> {
        let new_cell = self.transform_cell(&magnetic_cell.cell);
        MagneticCell::new(
            new_cell.lattice,
            new_cell.positions,
            new_cell.numbers,
            self.transform_magnetic_moments(&magnetic_cell.magnetic_moments),
        )
    }
}

impl Mul for UnimodularTransformation {
    type Output = Self;

    // (P_lhs, p_lhs) * (P_rhs, p_rhs) = (P_lhs * P_rhs, P_lhs * p_rhs + p_lhs)
    fn mul(self, rhs: Self) -> Self::Output {
        let new_linear = self.linear * rhs.linear;
        let new_origin_shift = self.linear.map(|e| e as f64) * rhs.origin_shift + self.origin_shift;
        Self::new(new_linear, new_origin_shift)
    }
}

/// Represent change of origin and basis for an affine space
#[derive(Debug, Clone)]
pub struct Transformation {
    pub linear: Linear,
    pub origin_shift: OriginShift,
    #[allow(dead_code)]
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

    #[allow(dead_code)]
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
    pub fn transform_operation(&self, operation: &Operation) -> Option<Operation> {
        transform_operation_as_f64(
            operation,
            &self.linear.map(|e| e as f64),
            &self.linear_inv,
            &self.origin_shift,
        )
    }

    /// (P, p)^-1 (W, w) (P, p)
    /// This function may decrease the number of operations if the transformation is not compatible with an operation.
    pub fn transform_operations(&self, operations: &[Operation]) -> Vec<Operation> {
        operations
            .iter()
            .filter_map(|ops| self.transform_operation(ops))
            .collect()
    }

    /// (P, p) (W, w) (P, p)^-1
    pub fn inverse_transform_operation(&self, operation: &Operation) -> Option<Operation> {
        transform_operation_as_f64(
            operation,
            &self.linear_inv,
            &self.linear.map(|e| e as f64),
            &(-self.linear_as_f64() * self.origin_shift),
        )
    }

    /// (P, p) (W, w) (P, p)^-1
    /// This function may decrease the number of operations if the transformation is not compatible with an operation.
    pub fn inverse_transform_operations(&self, operations: &[Operation]) -> Vec<Operation> {
        operations
            .iter()
            .filter_map(|ops| self.inverse_transform_operation(ops))
            .collect()
    }

    pub fn transform_magnetic_operation(
        &self,
        magnetic_operation: &MagneticOperation,
    ) -> Option<MagneticOperation> {
        let new_operation = self.transform_operation(&magnetic_operation.operation)?;
        Some(MagneticOperation::from_operation(
            new_operation,
            magnetic_operation.time_reversal,
        ))
    }

    pub fn transform_magnetic_operations(
        &self,
        magnetic_operations: &[MagneticOperation],
    ) -> Vec<MagneticOperation> {
        magnetic_operations
            .iter()
            .filter_map(|mops| self.transform_magnetic_operation(mops))
            .collect()
    }

    pub fn inverse_transform_magnetic_operation(
        &self,
        magnetic_operation: &MagneticOperation,
    ) -> Option<MagneticOperation> {
        let new_operation = self.inverse_transform_operation(&magnetic_operation.operation)?;
        Some(MagneticOperation::from_operation(
            new_operation,
            magnetic_operation.time_reversal,
        ))
    }

    pub fn inverse_transform_magnetic_operations(
        &self,
        magnetic_operations: &[MagneticOperation],
    ) -> Vec<MagneticOperation> {
        magnetic_operations
            .iter()
            .filter_map(|mops| self.inverse_transform_magnetic_operation(mops))
            .collect()
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

        let new_num_atoms = cell.num_atoms() * lattice_points.len();
        let mut new_positions = Vec::with_capacity(new_num_atoms);
        let mut new_numbers = Vec::with_capacity(new_num_atoms);
        let mut site_mapping = Vec::with_capacity(new_num_atoms);
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

    pub fn transform_magnetic_cell<M: MagneticMoment>(
        &self,
        magnetic_cell: &MagneticCell<M>,
    ) -> (MagneticCell<M>, Vec<usize>) {
        let (new_cell, site_mapping) = self.transform_cell(&magnetic_cell.cell);
        let new_magnetic_moments = site_mapping
            .iter()
            .map(|&i| magnetic_cell.magnetic_moments[i].clone()) // magnetic moments are not transformed
            .collect();
        (
            MagneticCell::new(
                new_cell.lattice,
                new_cell.positions,
                new_cell.numbers,
                new_magnetic_moments,
            ),
            site_mapping,
        )
    }
}

/// Transform operation (rotation, translation) by transformation (linear, origin_shift).
fn transform_operation_as_f64(
    operation: &Operation,
    linear: &Matrix3<f64>,
    linear_inv: &Matrix3<f64>,
    origin_shift: &OriginShift,
) -> Option<Operation> {
    let new_rotation =
        (linear_inv * operation.rotation.map(|e| e as f64) * linear).map(|e| e.round() as i32);

    // Check if `new_rotation` is an integer matrix
    let recovered =
        (linear * new_rotation.map(|e| e as f64) * linear_inv).map(|e| e.round() as i32);
    if recovered != operation.rotation {
        return None;
    }

    let new_translation = linear_inv
        * (operation.rotation.map(|e| e as f64) * origin_shift + operation.translation
            - origin_shift);
    Some(Operation::new(new_rotation, new_translation))
}

#[cfg(test)]
mod tests {
    use nalgebra::matrix;

    use super::Transformation;
    use crate::base::operation::{Operation, Translation};

    #[test]
    fn test_incompatible_transformation() {
        let transformation = Transformation::from_linear(matrix![
            1, 0, 0;
            0, 1, 0;
            0, 0, 2;
        ]);
        // threefold rotation
        let operation = Operation::new(
            matrix![
                0, 0, 1;
                1, 0, 0;
                0, 1, 0;
            ],
            Translation::zeros(),
        );
        assert!(transformation.transform_operation(&operation).is_none());
    }
}
