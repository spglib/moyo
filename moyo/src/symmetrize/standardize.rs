use log::warn;
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{Matrix3, Vector3, vector};
use std::collections::HashMap;

use super::conventional_cell::{
    AXIS_PERMUTATIONS3, monoclinic_candidate_corrections, monoclinic_rank_key,
    orthorhombic_rank_key, select_conventional_correction,
};
use super::wyckoff::{assign_wyckoffs_by_orbit, group_sites_by_orbit, match_wyckoff_coordinates};
use crate::base::{
    Cell, EPS, Lattice, Linear, MoyoError, Operations, Permutation, Position, Rotations,
    Transformation, UnimodularTransformation, project_rotations,
};
use crate::data::{
    HallNumber, HallSymbol, LatticeSystem, WyckoffPosition, arithmetic_crystal_class_entry,
    hall_symbol_entry, iter_wyckoff_positions,
};
use crate::identify::SpaceGroup;

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
            space_group,
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
        space_group: &SpaceGroup,
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
        let entry =
            hall_symbol_entry(space_group.hall_number).ok_or(MoyoError::StandardizationError)?;

        // Prepare operations in primitive standard
        let hs = HallSymbol::from_hall_number(space_group.hall_number)
            .ok_or(MoyoError::StandardizationError)?;
        let (conv_std_operations, prim_std_operations) = hs.traverse_and_primitive_traverse();

        // To standardized primitive cell
        let lattice_system = arithmetic_crystal_class_entry(entry.arithmetic_number)
            .unwrap()
            .lattice_system();
        // For monoclinic and orthorhombic systems, the identified setting leaves some
        // freedom in the conventional basis (the affine normalizer). The chosen
        // correction is folded into `prim_transformation` so that the primitive
        // standardized cell is always related to the conventional one by the fixed
        // centering matrix `entry.centering.linear()`.
        let conv_lattice_tmp = Transformation::from_linear(
            space_group.transformation.linear * entry.centering.linear(),
        )
        .transform_lattice(&prim_cell.lattice);
        let (prim_transformation, conv_trans_linear) = match lattice_system {
            LatticeSystem::Triclinic => (
                standardize_triclinic_cell(&prim_cell.lattice, &space_group.transformation),
                Linear::identity(),
            ),
            LatticeSystem::Monoclinic => {
                let prim_correction = select_conventional_correction(
                    &conv_lattice_tmp,
                    entry.centering,
                    &prim_std_operations,
                    &hs.primitive_generators(),
                    &monoclinic_candidate_corrections(&conv_lattice_tmp),
                    monoclinic_rank_key,
                    epsilon,
                );
                (
                    space_group.transformation.clone() * prim_correction,
                    entry.centering.linear(),
                )
            }
            LatticeSystem::Orthorhombic => {
                let prim_correction = select_conventional_correction(
                    &conv_lattice_tmp,
                    entry.centering,
                    &prim_std_operations,
                    &hs.primitive_generators(),
                    &AXIS_PERMUTATIONS3,
                    orthorhombic_rank_key,
                    epsilon,
                );
                (
                    space_group.transformation.clone() * prim_correction,
                    entry.centering.linear(),
                )
            }
            _ => (space_group.transformation.clone(), entry.centering.linear()),
        };

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

        // prim_transformation * (conv_trans_linear, 0)
        let transformation = Transformation::new(
            prim_transformation.linear * conv_trans_linear,
            prim_transformation.origin_shift,
        );
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

/// Niggli reduction for distorted triclinic lattice systems is numerically so challenging.
/// Thus, we skip checking reduction condition.
fn standardize_triclinic_cell(
    lattice: &Lattice,
    transformation_to_prim_std: &UnimodularTransformation,
) -> UnimodularTransformation {
    let lattice_prim_std_tmp = transformation_to_prim_std.transform_lattice(lattice);
    let (_, niggli_linear) = lattice_prim_std_tmp.unchecked_niggli_reduce();
    UnimodularTransformation::new(
        niggli_linear * transformation_to_prim_std.linear,
        transformation_to_prim_std.origin_shift,
    )
}

/// Symmetrize positions by site symmetry groups. Operates on a slice rather
/// than `&Cell` so the layer pipeline (which threads `LayerCell::positions()`)
/// can reuse it.
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

fn symmetrize_lattice(lattice: &Lattice, rotations: &Rotations) -> (Lattice, Matrix3<f64>) {
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
    use nalgebra::matrix;

    use super::symmetrize_lattice;
    use crate::base::{Lattice, traverse};
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
}
