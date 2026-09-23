use log::warn;
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{Matrix3, Vector3, vector};
use std::collections::HashMap;

use super::conventional_coordinate_system::ConventionalCoordinateSystem;
use super::wyckoff::{assign_wyckoffs_by_orbit, group_sites_by_orbit, match_wyckoff_coordinates};
use crate::base::{
    Cell, EPS, Lattice, MoyoError, Operations, Permutation, Position, Rotations, Transformation,
    UnimodularTransformation, project_rotations,
};
use crate::data::{HallNumber, WyckoffPosition, iter_wyckoff_positions};
use crate::identify::SpaceGroup;

/// Result of the full standardization pipeline for a primitive cell.
///
/// [`ConventionalCoordinateSystem`] selects the basis and origin. This type then
/// refines atomic positions, constructs the primitive and conventional cells,
/// applies the requested Cartesian orientation, and assigns Wyckoff positions.
pub struct StandardizedCell {
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    pub prim_cell: Cell,
    /// Transformation from the input primitive cell to the primitive standardized cell.
    pub prim_transformation: UnimodularTransformation,
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    pub cell: Cell,
    /// Wyckoff positions of sites in the `cell`
    pub wyckoffs: Vec<WyckoffPosition>,
    /// Transformation from the input primitive cell to the standardized cell.
    pub transformation: Transformation,
    /// Rotation matrix to map the lattice of the input primitive cell to that of the standardized cell.
    // ------------------------------------------------------------------------
    // Miscellaneous
    // ------------------------------------------------------------------------
    pub rotation_matrix: Matrix3<f64>,
    /// Mapping from the site in the `cell` to that in the `prim_cell`
    pub site_mapping: Vec<usize>,
}

impl StandardizedCell {
    /// Standardize the input **primitive** cell.
    /// For triclinic space groups, Niggli reduction is performed.
    /// Basis vectors are rotated to be a upper triangular matrix.
    pub fn new(
        prim_cell: &Cell,
        prim_operations: &Operations,
        prim_permutations: &[Permutation],
        space_group: &SpaceGroup,
        symprec: f64,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<Self, MoyoError> {
        let coordinate_system =
            ConventionalCoordinateSystem::new(&prim_cell.lattice, space_group, epsilon)?;
        let (
            prim_std_cell,
            prim_std_permutations,
            prim_transformation,
            std_cell,
            transformation,
            rotation_matrix,
            site_mapping,
        ) = Self::standardize_and_symmetrize_cell(
            prim_cell,
            prim_operations,
            prim_permutations,
            coordinate_system,
            epsilon,
            rotate_basis,
        )?;

        let wyckoffs = Self::assign_wyckoffs(
            &prim_std_cell,
            &prim_std_permutations,
            &std_cell,
            &site_mapping,
            space_group.hall_number,
            symprec,
        )?;

        Ok(StandardizedCell {
            // Primitive standardized cell
            prim_cell: prim_std_cell,
            prim_transformation,
            // Standardized cell
            cell: std_cell,
            wyckoffs,
            transformation,
            // Miscellaneous
            rotation_matrix,
            site_mapping,
        })
    }

    #[allow(clippy::type_complexity)]
    fn standardize_and_symmetrize_cell(
        prim_cell: &Cell,
        prim_operations: &Operations,
        prim_permutations: &[Permutation],
        coordinate_system: ConventionalCoordinateSystem,
        epsilon: f64,
        rotate_basis: bool,
    ) -> Result<
        (
            Cell,
            Vec<Permutation>,
            UnimodularTransformation,
            Cell,
            Transformation,
            Matrix3<f64>,
            Vec<usize>,
        ),
        MoyoError,
    > {
        let ConventionalCoordinateSystem {
            prim_transformation,
            conv_trans_linear,
            transformation,
            conv_std_operations,
            prim_std_operations,
        } = coordinate_system;

        let prim_std_cell_tmp = prim_transformation.transform_cell(prim_cell);

        // Symmetrize positions of prim_std_cell by refined symmetry operations.
        // Reorder permutations because prim_std_operations (from the Hall-symbol
        // traversal) is in a different order than `prim_operations` (from the
        // symmetry search).
        let prim_std_permutations = align_primitive_permutations(
            &prim_transformation,
            prim_operations,
            prim_permutations,
            &prim_std_operations,
        )?;
        let new_prim_std_positions = symmetrize_positions(
            &prim_std_cell_tmp.positions,
            &prim_std_operations,
            &prim_std_permutations,
            epsilon,
        );

        // Note: prim_transformation.transform_cell does not change the order of sites
        let prim_std_cell = Cell::new(
            prim_std_cell_tmp.lattice.clone(),
            new_prim_std_positions,
            prim_std_cell_tmp.numbers.clone(),
        );

        // To (conventional) standardized cell
        let (std_cell, site_mapping) =
            Transformation::from_linear(conv_trans_linear).transform_cell(&prim_std_cell);

        if rotate_basis {
            // Symmetrize lattice
            let (_, rotation_matrix) =
                symmetrize_lattice(&std_cell.lattice, &project_rotations(&conv_std_operations));
            Ok((
                prim_std_cell.rotate(&rotation_matrix),
                prim_std_permutations,
                prim_transformation.clone(),
                std_cell.rotate(&rotation_matrix),
                transformation,
                rotation_matrix,
                site_mapping,
            ))
        } else {
            Ok((
                prim_std_cell,
                prim_std_permutations,
                prim_transformation.clone(),
                std_cell,
                transformation,
                Matrix3::identity(),
                site_mapping,
            ))
        }
    }

    fn assign_wyckoffs(
        prim_std_cell: &Cell,
        prim_std_permutations: &[Permutation],
        std_cell: &Cell,
        site_mapping: &[usize],
        hall_number: HallNumber,
        symprec: f64,
    ) -> Result<Vec<WyckoffPosition>, MoyoError> {
        let group = group_sites_by_orbit(
            prim_std_cell.num_atoms(),
            prim_std_permutations,
            site_mapping,
            std_cell.num_atoms(),
        );
        assign_wyckoffs_by_orbit(&group, &std_cell.positions, |position, multiplicity| {
            iter_wyckoff_positions(hall_number, multiplicity)
                .find(|w| {
                    match_wyckoff_coordinates(position, w.coordinates, &std_cell.lattice, symprec)
                })
                .cloned()
        })
    }
}

/// Align a set of primitive-cell permutations (as produced by the symmetry
/// search) with a target operation list (as produced by traversing a Hall
/// symbol) by matching rotation matrices. Both pipelines need this because
/// the search and the database traversal generate operations in different
/// orders. Errors with `StandardizationError` if any target rotation is
/// missing from the input — i.e. the input set is not closed under the
/// transformation, which would indicate a primitive-cell or
/// hall-symbol-database bug rather than a user error.
pub(super) fn align_primitive_permutations(
    prim_transformation: &UnimodularTransformation,
    prim_operations: &Operations,
    prim_permutations: &[Permutation],
    target_operations: &Operations,
) -> Result<Vec<Permutation>, MoyoError> {
    let mut permutation_mapping = HashMap::new();
    let prim_rotations =
        project_rotations(&prim_transformation.transform_operations(prim_operations));
    for (rot, perm) in prim_rotations.iter().zip(prim_permutations.iter()) {
        permutation_mapping.insert(*rot, perm.clone());
    }
    target_operations
        .iter()
        .map(|ops| {
            permutation_mapping
                .get(&ops.rotation)
                .cloned()
                .ok_or(MoyoError::StandardizationError)
        })
        .collect()
}

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
