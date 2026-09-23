//! Refine lattice metrics and atomic positions under fixed symmetry operations.
//!
//! Callers supply the coordinate system, target operations, and aligned atom
//! permutations. These helpers return the refined quantities; the standardization
//! pipelines construct the output cells and choose their Cartesian orientation.

use log::warn;
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{Matrix3, Vector3, vector};

use crate::base::{EPS, Lattice, Operations, Permutation, Position, Rotations};

/// Refine fractional positions under fixed symmetry operations and atom mappings.
///
/// `operations` and `permutations` must have the same nonzero length, with
/// `permutations[k]` mapping the input sites under `operations[k]`. The operations
/// must form a group modulo lattice translations, and the mappings must respect
/// its action on the sites. Both operations and positions use the same basis and
/// origin. The caller supplies these correspondences; this function does not
/// search for symmetry or choose a coordinate system.
///
/// Average the mapped positions after reducing each fractional displacement
/// component by its nearest integer. The input must be close enough to the target
/// symmetry for these image choices to be consistent. Output positions retain
/// the input site order and are not wrapped into the unit cell. `epsilon` only controls
/// warnings about large fractional displacements; it does not stop refinement.
pub(super) fn symmetrize_positions(
    positions: &[Position],
    operations: &Operations,
    permutations: &[Permutation],
    epsilon: f64,
) -> Vec<Position> {
    // operations[k] maps site-i to site-permutations[k].apply(i)
    // Thus, it maps site-`inverse_permutations[k].apply(i)` to site-i.
    let inverse_permutations = permutations
        .iter()
        .map(|permutation| permutation.inverse())
        .collect::<Vec<_>>();

    (0..positions.len())
        .map(|i| {
            let mut acc = Vector3::zeros();
            for (inv_perm, operation) in inverse_permutations.iter().zip(operations.iter()) {
                let mut frac_displacements = operation.rotation.map(|e| e as f64)
                    * positions[inv_perm.apply(i)]
                    + operation.translation
                    - positions[i];
                frac_displacements -= frac_displacements.map(|e| e.round()); // in [-0.5, 0.5]
                acc += frac_displacements;
            }
            acc /= permutations.len() as f64;
            if acc.abs().max() > epsilon {
                warn!(
                    "Large displacement during symmetrization: {:?} for site {}",
                    acc, i
                )
            }
            positions[i] + acc
        })
        .collect::<Vec<_>>()
}

/// Refine a lattice metric under a fixed, nonempty crystallographic point group.
///
/// `rotations` must form a finite group expressed in the input lattice basis.
/// Average `W^T G W` over that group, where `G` is the input metric. Form an
/// upper-triangular basis from its Cholesky factor, then adjust the axis signs to
/// preserve the input handedness.
///
/// The second result is the proper orthogonal factor from the QR decomposition
/// of `A_refined A_input^-1`. It describes Cartesian orientation separately from
/// metric refinement: rotating the input by this matrix generally does not yield
/// the refined lattice. Fractional coordinates and the choice of basis and origin
/// are outside this function's responsibility.
pub(super) fn symmetrize_lattice(
    lattice: &Lattice,
    rotations: &Rotations,
) -> (Lattice, Matrix3<f64>) {
    let metric_tensor = lattice.metric_tensor();
    let mut symmetrized_metric_tensor: Matrix3<f64> = rotations
        .iter()
        .map(|rotation| {
            rotation.transpose().map(|e| e as f64) * metric_tensor * rotation.map(|e| e as f64)
        })
        .sum();
    symmetrized_metric_tensor /= rotations.len() as f64;

    // Upper-triangular basis
    let mut tri_basis = Cholesky::new_unchecked(symmetrized_metric_tensor)
        .l()
        .transpose();
    // Remove axis-direction freedom
    let diagonal_signs = Matrix3::<f64>::from_diagonal(&vector![
        sign(tri_basis[(0, 0)]),
        sign(tri_basis[(1, 1)]),
        sign(tri_basis[(2, 2)])
    ]);
    tri_basis *= diagonal_signs;
    // Adjust handedness
    if sign(lattice.basis.determinant()) * sign(tri_basis.determinant()) < 0.0 {
        tri_basis *= Matrix3::<f64>::from_diagonal(&vector![1.0, 1.0, -1.0]);
    }

    // tri_basis \approx orthogonal_matrix * lattice.basis
    // QR(tri_basis * lattice.basis^-1) = rotation_matrix * strain
    let mut rotation_matrix = QR::new(tri_basis * lattice.basis.try_inverse().unwrap()).q();
    if rotation_matrix.determinant() < 0.0 {
        rotation_matrix *= -1.0;
    }

    (Lattice::new(tri_basis.transpose()), rotation_matrix)
}

fn sign(x: f64) -> f64 {
    if x > EPS {
        1.0
    } else if x < -EPS {
        -1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::{Matrix3, Rotation3, Vector3, matrix, vector};

    use super::{symmetrize_lattice, symmetrize_positions};
    use crate::base::{Lattice, Operation, Permutation, traverse};
    use crate::data::{GeometricCrystalClass, PointGroupRepresentative};

    #[test]
    fn test_symmetrize_lattice_cubic() {
        let lattice = Lattice::new(matrix![
            1.0, 0.0, 0.0001;
            0.0, -0.999, 0.0;
            0.0, 0.0, -1.0001;
        ]);
        let rep = PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::Oh);
        let rotations = traverse(&rep.generators);

        let (new_lattice, rotation_matrix) = symmetrize_lattice(&lattice, &rotations);
        assert_relative_eq!(new_lattice.basis[(1, 1)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(2, 2)], new_lattice.basis[(0, 0)]);
        assert_relative_eq!(new_lattice.basis[(0, 1)], 0.0);
        assert_relative_eq!(new_lattice.basis[(0, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(1, 2)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 0)], 0.0);
        assert_relative_eq!(new_lattice.basis[(2, 1)], 0.0);

        assert_relative_eq!(
            rotation_matrix * lattice.basis,
            new_lattice.basis,
            epsilon = 1e-2
        );
    }

    #[rstest::rstest]
    #[case::right_handed(1.0)]
    #[case::left_handed(-1.0)]
    fn test_symmetrize_lattice_metric_and_orientation(#[case] handedness: f64) {
        let lattice = Lattice::new(matrix![
            1.001, 0.002, 0.0;
            0.0, 0.999, -0.001;
            0.0, 0.0, handedness * 1.0001;
        ])
        .rotate(Rotation3::from_euler_angles(0.2, -0.3, 0.4).matrix());
        let rep = PointGroupRepresentative::from_geometric_crystal_class(GeometricCrystalClass::Oh);
        let rotations = traverse(&rep.generators);

        let (refined, rotation) = symmetrize_lattice(&lattice, &rotations);

        // Cubic averaging preserves the trace and makes all three lengths equal.
        let squared_length = lattice.metric_tensor().trace() / 3.0;
        let length = squared_length.sqrt();
        assert_relative_eq!(
            refined.basis,
            Matrix3::from_diagonal(&vector![length, length, handedness * length]),
            epsilon = 1e-12
        );
        assert!((refined.metric_tensor() - lattice.metric_tensor()).norm() > 1e-3);
        for rotation in &rotations {
            let rotation = rotation.map(|e| e as f64);
            assert_relative_eq!(
                rotation.transpose() * refined.metric_tensor() * rotation,
                refined.metric_tensor(),
                epsilon = 1e-12
            );
        }
        assert_relative_eq!(
            rotation.transpose() * rotation,
            Matrix3::identity(),
            epsilon = 1e-12
        );
        assert_relative_eq!(rotation.determinant(), 1.0, epsilon = 1e-12);

        // The refined lattice is a fixed point of the metric projection.
        let (refined_again, _) = symmetrize_lattice(&refined, &rotations);
        assert_relative_eq!(refined_again.basis, refined.basis, epsilon = 1e-12);
    }

    #[test]
    fn test_symmetrize_positions_screw_orbit() {
        let screw = Operation::new(
            matrix![0, -1, 0; 1, -1, 0; 0, 0, 1],
            vector![0.0, 0.0, 1.0 / 3.0],
        );
        let operations = vec![screw.clone() * screw.clone(), Operation::identity(), screw];
        let permutations = vec![
            Permutation::new(vec![2, 0, 1]),
            Permutation::identity(3),
            Permutation::new(vec![1, 2, 0]),
        ];
        // Perturb a screw orbit, keeping one site outside the unit cell.
        let positions = vec![
            vector![1.2, 0.3, 0.95] + vector![0.001, -0.002, 0.003],
            vector![0.7, 0.9, 0.95 + 1.0 / 3.0 - 1.0] + vector![-0.004, 0.005, -0.006],
            vector![0.1, 0.8, 0.95 + 2.0 / 3.0 - 1.0] + vector![0.007, -0.008, 0.009],
        ];
        let expected = [
            vector![1.206, 0.3 + 0.017 / 3.0, 0.952],
            vector![
                0.7 - 0.017 / 3.0,
                0.9 + 0.001 / 3.0,
                0.952 + 1.0 / 3.0 - 1.0
            ],
            vector![0.1 - 0.001 / 3.0, 0.794, 0.952 + 2.0 / 3.0 - 1.0],
        ];

        let refined = symmetrize_positions(&positions, &operations, &permutations, 1e-2);
        assert_eq!(refined.len(), positions.len());
        for (position, expected) in refined.iter().zip(expected.iter()) {
            assert_relative_eq!(position, expected, epsilon = 1e-12);
        }
        for (operation, permutation) in operations.iter().zip(permutations.iter()) {
            for (i, position) in refined.iter().enumerate() {
                let displacement = operation.rotation.map(|e| e as f64) * position
                    + operation.translation
                    - refined[permutation.apply(i)];
                assert_relative_eq!(
                    displacement - displacement.map(|e| e.round()),
                    Vector3::zeros(),
                    epsilon = 1e-12
                );
            }
        }
        let refined_again = symmetrize_positions(&refined, &operations, &permutations, 1e-2);
        for (position, refined_again) in refined.iter().zip(refined_again.iter()) {
            assert_relative_eq!(position, refined_again, epsilon = 1e-12);
        }
    }
}
