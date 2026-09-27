use std::ops::{Deref, Mul};

use itertools::iproduct;
use nalgebra::base::{Matrix3, Vector3};
use serde::Serialize;

use super::cell::Cell;
use super::lattice::Lattice;
use super::layer::{LayerCell, LayerLattice};
use super::magnetic_cell::{MagneticCell, MagneticMoment};
use super::operation::{MagneticOperation, Operation};
use super::tolerance::EPS;
use crate::math::SNF;
use crate::utils::{to_3_slice, to_3x3_slice};

pub type UnimodularLinear = Matrix3<i32>;
pub type Linear = Matrix3<i32>;
/// Origin shift in a crystallographic basis
pub type OriginShift = Vector3<f64>;

/// Representatives of `Z^3 / linear * Z^3` in the input basis.
pub(crate) fn lattice_points(linear: &Linear) -> Vec<Vector3<i32>> {
    // Consider the Smith normal form D = L * linear * R. Two lattice points
    // n and n' are equivalent exactly when L * n = L * n' (mod D * Z^3).
    // Thus L^-1 maps the product of the diagonal residue classes to one
    // representative of every coset.
    let snf = SNF::new(linear);
    let linear_inverse = snf
        .l
        .map(|element| element as f64)
        .try_inverse()
        .expect("Smith transformation is unimodular")
        .map(|element| element.round() as i32);

    iproduct!(0..snf.d[(0, 0)], 0..snf.d[(1, 1)], 0..snf.d[(2, 2)])
        .map(|(f0, f1, f2)| linear_inverse * Vector3::new(f0, f1, f2))
        .collect()
}

/// Change of origin and primitive basis with determinant +1 or -1.
///
/// The linear part and origin shift are immutable. The linear part's immutability
/// preserves unimodularity and its cached inverse.
///
/// ```compile_fail
/// use moyo::base::UnimodularTransformation;
/// use nalgebra::Matrix3;
/// let mut transformation = UnimodularTransformation::from_linear(Matrix3::identity());
/// transformation.linear = Matrix3::zeros();
/// ```
///
/// ```compile_fail
/// use moyo::base::UnimodularTransformation;
/// use nalgebra::Vector3;
/// let mut transformation = UnimodularTransformation::from_origin_shift(Vector3::zeros());
/// transformation.origin_shift = Vector3::zeros();
/// ```
#[derive(Debug, Clone, Serialize)]
pub struct UnimodularTransformation {
    linear: UnimodularLinear,
    origin_shift: OriginShift,
    // Inverse of unimodular matrix is also unimodular
    linear_inv: UnimodularLinear,
}

impl UnimodularTransformation {
    /// Construct a transformation, panicking unless the determinant is +/- 1
    /// and the integer inverse fits in `i32`.
    pub fn new(linear: UnimodularLinear, origin_shift: OriginShift) -> Self {
        // i128 accommodates products of three i32 entries without rounding.
        let wide = linear.map(i128::from);
        let cofactors = Matrix3::from_columns(&[
            wide.column(1).cross(&wide.column(2)),
            wide.column(2).cross(&wide.column(0)),
            wide.column(0).cross(&wide.column(1)),
        ]);
        let det = wide.column(0).dot(&cofactors.column(0));
        if det.abs() != 1 {
            panic!("Determinant of unimodular transformation must be +/- 1.");
        }

        let linear_inv = cofactors.transpose().map(|e| {
            i32::try_from(e / det).expect("Inverse of unimodular transformation must fit in i32.")
        });

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

    /// Immutable linear part of this transformation.
    pub fn linear(&self) -> &UnimodularLinear {
        &self.linear
    }

    /// Immutable origin shift of this transformation.
    pub fn origin_shift(&self) -> &OriginShift {
        &self.origin_shift
    }

    /// Exact determinant of the linear part, either +1 or -1.
    pub fn determinant(&self) -> i32 {
        let wide = self.linear.map(i128::from);
        wide.column(0).dot(&wide.column(1).cross(&wide.column(2))) as i32
    }

    pub fn linear_as_f64(&self) -> Matrix3<f64> {
        self.linear.map(|e| e as f64)
    }

    /// Returns the linear part as a 3x3 integer array.
    pub fn linear_as_array(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.linear)
    }

    /// Returns the origin shift as a `[f64; 3]` array.
    pub fn origin_shift_as_array(&self) -> [f64; 3] {
        to_3_slice(&self.origin_shift)
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

    /// Apply this transformation to a layer cell. Layer-shape preservation is
    /// the caller's responsibility: this is `pub(crate)` so only in-crate
    /// helpers, which know the linear part has the layer block form
    /// `W_i3 = W_3i = 0`, `|W_33| = 1`, can call it.
    pub(crate) fn transform_layer_cell(&self, cell: &LayerCell) -> LayerCell {
        let new_bulk = self.transform_cell(&cell.as_cell());
        LayerCell::new_unchecked(
            LayerLattice::new_unchecked(new_bulk.lattice),
            new_bulk.positions,
            new_bulk.numbers,
        )
    }

    pub fn transform_magnetic_moments<M: MagneticMoment>(&self, magnetic_moments: &[M]) -> Vec<M> {
        // Magnetic moments are not transformed
        magnetic_moments.to_owned()
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

/// Orientation-preserving change of origin and primitive basis (determinant +1).
///
/// Coordinate transformations are shared with [`UnimodularTransformation`] through
/// immutable dereferencing. Inversion and composition of proper transformations
/// return this type; conversion to the general type is infallible.
///
/// ```compile_fail
/// use moyo::base::ProperUnimodularTransformation;
/// use nalgebra::Matrix3;
/// let mut transformation = ProperUnimodularTransformation::from_linear(Matrix3::identity());
/// transformation.linear()[(0, 0)] = -1;
/// ```
#[derive(Debug, Clone, Serialize)]
#[serde(transparent)]
pub struct ProperUnimodularTransformation(UnimodularTransformation);

/// An orientation-reversing transformation cannot be converted to a proper one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("Determinant of proper unimodular transformation must be +1.")]
pub struct ImproperTransformationError;

impl ProperUnimodularTransformation {
    /// Construct a transformation, panicking unless the determinant is +1
    /// and the integer inverse fits in `i32`.
    pub fn new(linear: UnimodularLinear, origin_shift: OriginShift) -> Self {
        Self::try_from(UnimodularTransformation::new(linear, origin_shift))
            .expect("Determinant of proper unimodular transformation must be +1.")
    }

    pub fn from_linear(linear: UnimodularLinear) -> Self {
        Self::new(linear, OriginShift::zeros())
    }

    pub fn from_origin_shift(origin_shift: OriginShift) -> Self {
        Self::new(UnimodularLinear::identity(), origin_shift)
    }

    pub fn inverse(&self) -> Self {
        Self(self.0.inverse())
    }
}

impl Deref for ProperUnimodularTransformation {
    type Target = UnimodularTransformation;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<ProperUnimodularTransformation> for UnimodularTransformation {
    fn from(value: ProperUnimodularTransformation) -> Self {
        value.0
    }
}

impl TryFrom<UnimodularTransformation> for ProperUnimodularTransformation {
    type Error = ImproperTransformationError;

    fn try_from(value: UnimodularTransformation) -> Result<Self, Self::Error> {
        if value.determinant() == 1 {
            Ok(Self(value))
        } else {
            Err(ImproperTransformationError)
        }
    }
}

impl Mul for ProperUnimodularTransformation {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self(self.0 * rhs.0)
    }
}

impl Mul<ProperUnimodularTransformation> for UnimodularTransformation {
    type Output = Self;

    fn mul(self, rhs: ProperUnimodularTransformation) -> Self::Output {
        self * rhs.0
    }
}

impl Mul<UnimodularTransformation> for ProperUnimodularTransformation {
    type Output = UnimodularTransformation;

    fn mul(self, rhs: UnimodularTransformation) -> Self::Output {
        self.0 * rhs
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

    /// Returns the linear part as a 3x3 integer array.
    pub fn linear_as_array(&self) -> [[i32; 3]; 3] {
        to_3x3_slice(&self.linear)
    }

    /// Returns the origin shift as a `[f64; 3]` array.
    pub fn origin_shift_as_array(&self) -> [f64; 3] {
        to_3_slice(&self.origin_shift)
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

        let lattice_points = lattice_points(&self.linear);

        let new_num_atoms = cell.num_atoms() * lattice_points.len();
        let mut new_positions = Vec::with_capacity(new_num_atoms);
        let mut new_numbers = Vec::with_capacity(new_num_atoms);
        let mut site_mapping = Vec::with_capacity(new_num_atoms);
        for (i, (pos, number)) in cell.positions.iter().zip(cell.numbers.iter()).enumerate() {
            for lattice_point in lattice_points.iter() {
                // Fractional coordinates in the new sublattice: P^-1 (x + n - p).
                let new_position = (self.linear_inv
                    * (pos + lattice_point.map(|element| element as f64) - self.origin_shift))
                    .map(|e| e % 1.);
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

    /// Apply this transformation to a layer cell. The layer block form of
    /// `linear` (`W_i3 = W_3i = 0`, `|W_33| = 1`) is the caller's
    /// responsibility, so this is `pub(crate)` and meant for the layer
    /// pipeline's centering / standardization helpers only.
    pub(crate) fn transform_layer_cell(&self, cell: &LayerCell) -> (LayerCell, Vec<usize>) {
        let (new_bulk, site_mapping) = self.transform_cell(&cell.as_cell());

        // Bulk `transform_cell` reduces every fractional component by `% 1.`,
        // which would clip aperiodic stacking heights. Restore each
        // transformed `z` from the source site: the layer block form
        // (`W_i3 = W_3i = 0`, `|W_33| = 1`) forces `D_2 = 1` in the SNF,
        // so the sublattice points have `z = 0` and the z-projection is
        // `new_z = (1/W_33) * old_z = +/- old_z`.
        let w33_inv = 1.0 / (self.linear[(2, 2)] as f64);
        let mut new_positions = new_bulk.positions;
        for (new_pos, &orig_idx) in new_positions.iter_mut().zip(site_mapping.iter()) {
            new_pos[2] = w33_inv * cell.positions()[orig_idx][2];
        }

        let new_layer = LayerCell::new_unchecked(
            LayerLattice::new_unchecked(new_bulk.lattice),
            new_positions,
            new_bulk.numbers,
        );
        (new_layer, site_mapping)
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
    let rotation = linear_inv * operation.rotation.map(|e| e as f64) * linear;
    let new_rotation = rotation.map(|e| e.round() as i32);

    // Test integrality in the target basis: conjugating the rounding error back
    // to the source basis can hide a non-integer rotation.
    if (rotation - new_rotation.map(|e| e as f64)).abs().max() > EPS {
        return None;
    }

    let new_translation = linear_inv
        * (operation.rotation.map(|e| e as f64) * origin_shift + operation.translation
            - origin_shift);
    Some(Operation::new(new_rotation, new_translation))
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector};

    use super::{
        ImproperTransformationError, Lattice, LayerCell, ProperUnimodularTransformation,
        Transformation, UnimodularTransformation,
    };
    use crate::base::AngleTolerance;
    use crate::base::cell::Cell;
    use crate::base::operation::{Operation, Translation};

    #[test]
    fn test_proper_affine_composition_and_inverse() {
        let left = ProperUnimodularTransformation::new(
            matrix![1, 1, 0; 0, 1, 0; 0, 0, 1],
            vector![0.25, 0.5, 0.0],
        );
        let right = ProperUnimodularTransformation::new(
            matrix![1, 0, 0; 1, 1, 0; 0, 0, 1],
            vector![0.5, 0.0, 0.25],
        );
        let product: ProperUnimodularTransformation = left.clone() * right.clone();
        assert_eq!(*product.linear(), matrix![2, 1, 0; 1, 1, 0; 0, 0, 1]);
        assert_relative_eq!(*product.origin_shift(), vector![0.75, 0.5, 0.25]);
        assert_eq!(product.determinant(), 1);
        let inverse: ProperUnimodularTransformation = product.inverse();
        assert_eq!(*inverse.linear(), matrix![1, -1, 0; -1, 2, 0; 0, 0, 1]);
        assert_relative_eq!(*inverse.origin_shift(), vector![-0.25, -0.25, -0.25]);
        let cell = Cell::new(
            Lattice::new(matrix![3.0, 0.0, 0.0; 0.2, 4.0, 0.0; 0.1, 0.3, 5.0]),
            vec![vector![0.13, 0.27, 0.41]],
            vec![14],
        );
        let sequential = right.transform_cell(&left.transform_cell(&cell));
        let transformed = product.transform_cell(&cell);
        assert_relative_eq!(transformed.lattice.basis, sequential.lattice.basis);
        assert_relative_eq!(transformed.positions[0], sequential.positions[0]);
        let back = inverse.transform_cell(&transformed);
        assert_relative_eq!(back.lattice.basis, cell.lattice.basis);
        assert_relative_eq!(back.positions[0], cell.positions[0]);
        assert_eq!(back.numbers, cell.numbers);

        let operation =
            Operation::new(matrix![0, -1, 0; 1, 0, 0; 0, 0, 1], vector![0.0, 0.0, 0.25]);
        let sequential = right.transform_operation(&left.transform_operation(&operation));
        let transformed = product.transform_operation(&operation);
        assert_eq!(transformed.rotation, sequential.rotation);
        assert_relative_eq!(transformed.translation, sequential.translation);
    }

    #[test]
    fn test_unimodular_type_conversions_and_mixed_composition() {
        let proper = ProperUnimodularTransformation::from_origin_shift(vector![0.25, 0.0, 0.0]);
        let general = UnimodularTransformation::from(proper.clone());
        assert_eq!(
            serde_json::to_value(&proper).unwrap(),
            serde_json::to_value(&general).unwrap()
        );
        let restored = ProperUnimodularTransformation::try_from(general).unwrap();
        assert_eq!(restored.linear(), proper.linear());
        assert_relative_eq!(*restored.origin_shift(), *proper.origin_shift());

        let mirror = UnimodularTransformation::from_linear(matrix![-1, 0, 0; 0, 1, 0; 0, 0, 1]);
        assert_eq!(
            ProperUnimodularTransformation::try_from(mirror.clone()).unwrap_err(),
            ImproperTransformationError
        );
        let left: UnimodularTransformation = mirror.clone() * proper.clone();
        let right: UnimodularTransformation = proper * mirror.clone();
        assert_eq!(left.determinant(), -1);
        assert_eq!(right.determinant(), -1);
        assert_relative_eq!(*left.origin_shift(), vector![-0.25, 0.0, 0.0]);
        assert_relative_eq!(*right.origin_shift(), vector![0.25, 0.0, 0.0]);
        let even: UnimodularTransformation = mirror.clone() * mirror;
        assert!(ProperUnimodularTransformation::try_from(even).is_ok());
    }

    #[rstest::rstest]
    #[case(matrix![-1, 0, 0; 0, 1, 0; 0, 0, 1])]
    #[case(matrix![2, 0, 0; 0, 1, 0; 0, 0, 1])]
    #[case(matrix![0, 0, 0; 0, 1, 0; 0, 0, 1])]
    #[should_panic]
    fn test_proper_transformation_rejects_invalid_determinant(
        #[case] linear: nalgebra::Matrix3<i32>,
    ) {
        ProperUnimodularTransformation::from_linear(linear);
    }

    #[rstest::rstest]
    #[case(1)]
    #[case(-1)]
    fn test_unimodular_inverse_is_exact(#[case] sign: i32) {
        // The two O(n^2) terms in the determinant differ by one, below f64 precision.
        let n = 100_000_000;
        let linear = matrix![sign * n, n - 1, 0; sign * (n + 1), n, 0; 0, 0, 1];
        let transformation = UnimodularTransformation::from_linear(linear);
        assert_eq!(transformation.determinant(), sign);
        assert_eq!(
            *transformation.inverse().linear(),
            matrix![sign * n, sign * (1 - n), 0; -n - 1, n, 0; 0, 0, 1]
        );
        assert_eq!(*transformation.inverse().inverse().linear(), linear);
    }

    #[test]
    #[should_panic(expected = "Inverse of unimodular transformation must fit in i32.")]
    fn test_unimodular_transformation_rejects_inverse_overflow() {
        UnimodularTransformation::from_linear(matrix![1, 50_000, 0; 0, 1, 50_000; 0, 0, 1]);
    }

    #[test]
    fn test_unimodular_transformation_accepts_det_minus_one() {
        // Mirror in x: det = -1. Must construct without panicking.
        let mirror = matrix![
            -1, 0, 0;
             0, 1, 0;
             0, 0, 1;
        ];
        let t = UnimodularTransformation::from_linear(mirror);
        assert_eq!(t.linear_as_f64().determinant().round() as i32, -1);
        // Inverse of an involution is itself.
        let inv = t.inverse();
        assert_eq!(inv.linear, mirror);
        // Conjugating a 4-fold about z by the mirror yields a 4^-1 about z
        // (the rotation sense is reversed). Either way, the operation must
        // be a valid integer matrix; just check round-trip via inverse.
        let op = Operation::new(
            matrix![
                0, -1, 0;
                1,  0, 0;
                0,  0, 1;
            ],
            vector![0.25, 0.0, 0.0],
        );
        let conjugated = t.transform_operation(&op);
        let back = inv.transform_operation(&conjugated);
        assert_eq!(back.rotation, op.rotation);
        assert_relative_eq!(back.translation, op.translation);
    }

    #[test]
    #[should_panic(expected = "Determinant of unimodular transformation must be +/- 1.")]
    fn test_unimodular_transformation_rejects_det_two() {
        let _ = UnimodularTransformation::from_linear(matrix![
            2, 0, 0;
            0, 1, 0;
            0, 0, 1;
        ]);
    }

    #[test]
    #[should_panic(expected = "Determinant of unimodular transformation must be +/- 1.")]
    fn test_unimodular_transformation_rejects_singular() {
        let _ = UnimodularTransformation::from_linear(matrix![
            1, 0, 0;
            0, 1, 0;
            0, 0, 0;
        ]);
    }

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

    #[test]
    fn test_incompatible_transformation_with_integer_round_trip() {
        let transformation = Transformation::from_linear(matrix![
            1, 0, 0;
            0, 3, 1;
            0, 0, 2;
        ]);
        let operation = Operation::new(
            matrix![
                -1, 0, 0;
                 0, 0, 1;
                 0, -1, 0;
            ],
            Translation::zeros(),
        );
        // The conjugate contains 1/2 and 5/6. Rounding it and conjugating back
        // nevertheless recovers the original integer rotation after rounding.
        assert!(transformation.transform_operation(&operation).is_none());
    }

    #[test]
    fn test_transform_cell_with_origin_shift() {
        let transformation = Transformation::new(
            matrix![2, 1, 0; 0, 1, 0; 0, 0, 1],
            vector![0.37, -0.23, 1.19],
        );
        let cell = Cell::new(
            Lattice::new(matrix![3.0, 0.0, 0.0; 0.2, 4.0, 0.0; 0.1, 0.3, 5.0]),
            vec![vector![1.9, -0.96, 1.6], vector![-1.0, 1.94, -1.14]],
            vec![14, 8],
        );
        let (transformed, mapping) = transformation.transform_cell(&cell);
        assert_eq!(transformed.num_atoms(), 4);
        assert_relative_eq!(
            transformed.lattice.basis,
            cell.lattice.basis * transformation.linear_as_f64(),
            epsilon = 1e-12
        );
        // The doubled, sheared cell contains two images of each species.
        // These coordinates include the origin shift before changing basis.
        for (number, expected) in [
            (14, vector![0.13, 0.27, 0.41]),
            (14, vector![0.63, 0.27, 0.41]),
            (8, vector![0.23, 0.17, 0.67]),
            (8, vector![0.73, 0.17, 0.67]),
        ] {
            let matches = transformed
                .positions
                .iter()
                .enumerate()
                .filter(|(j, position)| {
                    let delta = *position - expected;
                    transformed.numbers[*j] == number
                        && (delta - delta.map(f64::round)).norm() < 1e-12
                })
                .collect::<Vec<_>>();
            assert_eq!(
                matches.len(),
                1,
                "missing image of species {number}: {expected:?}"
            );
            for (j, _) in matches {
                assert_eq!(cell.numbers[mapping[j]], number);
            }
        }
    }

    #[test]
    fn test_transform_layer_cell_preserves_aperiodic_z() {
        // Centering transform with layer block form (`W_33 = 1`,
        // `W_i3 = W_3i = 0`). The third axis is aperiodic for layer cells, so
        // z must be preserved verbatim even when it lies outside `[0, 1)`.
        let centering = matrix![
            1, -1, 0;
            1,  1, 0;
            0,  0, 1;
        ];
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]);
        // z = 1.7 is outside `[0, 1)` -- a thicker slab. Bulk transform's
        // `% 1.` would fold it to 0.7.
        let z_outside = 1.7;
        let cell = Cell::new(lattice, vec![vector![0.1, 0.2, z_outside]], vec![1]);
        let layer_cell = LayerCell::new(cell, 1e-4, AngleTolerance::Default).unwrap();

        let (transformed, _) =
            Transformation::from_linear(centering).transform_layer_cell(&layer_cell);

        // Centered -> primitive doubles the cell, so each input atom maps to
        // two output sites; both must keep `z = 1.7`.
        for pos in transformed.positions() {
            assert!((pos[2] - z_outside).abs() < 1e-12, "z wrapped: {}", pos[2]);
        }
    }
}
