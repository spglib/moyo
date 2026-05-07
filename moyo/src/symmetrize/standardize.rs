use itertools::{Itertools, iproduct};
use log::{debug, warn};
use nalgebra::linalg::{Cholesky, QR};
use nalgebra::{Matrix3, Vector3, vector};
use once_cell::sync::Lazy;
use std::collections::HashMap;

use crate::base::{
    Cell, EPS, Lattice, Linear, MoyoError, Operations, Permutation, Position, Rotations,
    Transformation, UnimodularTransformation, orbits_from_permutations, project_rotations,
};
use crate::data::{
    Centering, HallNumber, HallSymbol, LatticeSystem, WyckoffPosition, WyckoffPositionSpace,
    arithmetic_crystal_class_entry, hall_symbol_entry, iter_wyckoff_positions,
};
use crate::identify::SpaceGroup;
use crate::math::SNF;

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
        let (prim_transformation, conv_trans_linear) = match lattice_system {
            LatticeSystem::Triclinic => (
                standardize_triclinic_cell(&prim_cell.lattice, &space_group.transformation),
                Linear::identity(),
            ),
            LatticeSystem::Monoclinic => {
                let trans_std_prim_to_conv = standardize_monoclinic_conv_cell(
                    &prim_cell.lattice,
                    &space_group.transformation,
                    entry.centering,
                    &hs.generators,
                    epsilon,
                );
                (space_group.transformation.clone(), trans_std_prim_to_conv)
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

/// * `prim_num_atoms` - Number of atoms in the primitive cell
/// * `prim_permutations` - Permutations of the primitive cell
/// * `site_mapping` - Mapping from the site in cell to that in its primitive cell
pub fn orbits_in_cell(
    prim_num_atoms: usize,
    prim_permutations: &[Permutation],
    site_mapping: &[usize],
) -> Vec<usize> {
    // prim_orbits: [prim_num_atoms] -> [prim_num_atoms]
    let prim_orbits = orbits_from_permutations(prim_num_atoms, prim_permutations);

    let num_atoms = site_mapping.len();
    let mut map = HashMap::new();
    let mut orbits = vec![]; // [num_atoms] -> [num_atoms]
    for i in 0..num_atoms {
        // site_mapping: [num_atoms] -> [prim_num_atoms]
        let key = prim_orbits[site_mapping[i]]; // in [prim_num_atoms]
        map.entry(key).or_insert(i);
        orbits.push(*map.get(&key).unwrap());
    }
    orbits
}

/// Group sites of a conventional cell by crystallographic orbit. Shared
/// between the bulk and layer pipelines.
pub(super) struct OrbitGrouping {
    /// `[std_num_atoms] -> [num_orbits]`: orbit index for each site.
    pub mapping: Vec<usize>,
    /// `[num_orbits] -> [std_num_atoms]`: representative site for each orbit
    /// (the lowest-index site in the orbit). Used only for diagnostics.
    pub remapping: Vec<usize>,
    /// Number of sites per orbit, indexed by orbit. `multiplicities.len()` is
    /// the orbit count.
    pub multiplicities: Vec<usize>,
}

pub(super) fn group_sites_by_orbit(
    prim_num_atoms: usize,
    prim_permutations: &[Permutation],
    site_mapping: &[usize],
    std_num_atoms: usize,
) -> OrbitGrouping {
    let orbits = orbits_in_cell(prim_num_atoms, prim_permutations, site_mapping);

    let mut num_orbits = 0;
    let mut mapping = vec![0; std_num_atoms];
    let mut remapping = vec![];
    for i in 0..std_num_atoms {
        if orbits[i] == i {
            mapping[i] = num_orbits;
            remapping.push(i);
            num_orbits += 1;
        } else {
            mapping[i] = mapping[orbits[i]];
        }
    }

    let mut multiplicities = vec![0; num_orbits];
    for i in 0..std_num_atoms {
        multiplicities[mapping[i]] += 1;
    }

    OrbitGrouping {
        mapping,
        remapping,
        multiplicities,
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

/// Run the per-orbit Wyckoff-assignment loop over a precomputed
/// `OrbitGrouping`. The closure produces the database-specific Wyckoff
/// representative for an `(orbit_position, multiplicity)` pair (returning
/// `None` if no orbit in the database matches that multiplicity at that
/// position). Shared between bulk and layer pipelines: the only thing that
/// differs is the database lookup, which the caller passes in as the
/// closure.
pub(super) fn assign_wyckoffs_by_orbit<W, F>(
    group: &OrbitGrouping,
    positions: &[Position],
    mut find_for_orbit: F,
) -> Result<Vec<W>, MoyoError>
where
    W: Clone,
    F: FnMut(&Position, usize) -> Option<W>,
{
    let mut representative_wyckoffs: Vec<Option<W>> = vec![None; group.multiplicities.len()];
    for (i, position) in positions.iter().enumerate() {
        let orbit = group.mapping[i];
        if representative_wyckoffs[orbit].is_some() {
            continue;
        }
        if let Some(w) = find_for_orbit(position, group.multiplicities[orbit]) {
            representative_wyckoffs[orbit] = Some(w);
        }
    }

    for (orbit, wyckoff) in representative_wyckoffs.iter().enumerate() {
        if wyckoff.is_none() {
            debug!(
                "Failed to assign Wyckoff position with multiplicity {} at representative site {}",
                group.multiplicities[orbit], group.remapping[orbit]
            );
        }
    }
    let representative_wyckoffs = representative_wyckoffs
        .into_iter()
        .map(|w| w.ok_or(MoyoError::WyckoffPositionAssignmentError))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(group
        .mapping
        .iter()
        .map(|&orbit| representative_wyckoffs[orbit].clone())
        .collect())
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

/// Candidate unimodular corrections with entries in `[-1, 1]`.
///
/// This bounded search is used as a best-effort search space for monoclinic
/// conventional-cell standardization. Exhaustiveness is not claimed here.
static UNIMODULAR3_RANGE1: Lazy<Vec<UnimodularTransformation>> = Lazy::new(|| {
    (0..9)
        .map(|_| -1_i32..=1_i32)
        .multi_cartesian_product()
        .filter_map(|v| {
            let mat = Matrix3::new(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
            let det = mat.map(|e| e as f64).determinant().round() as i32;
            if det == 1 {
                Some(UnimodularTransformation::from_linear(mat))
            } else {
                None
            }
        })
        .collect()
});

/// Standardize monoclinic cell by choosing reduced basis vectors for perpendicular plane to the unique axis while keeping matrix representations of `conv_std_generators`.
fn standardize_monoclinic_conv_cell(
    prim_lattice: &Lattice,
    transformation_to_prim_std: &UnimodularTransformation,
    centering: Centering,
    conv_std_generators: &Operations,
    epsilon: f64,
) -> Linear {
    let mut candidate_conv_transformations = vec![];
    // Search candidate corrections in the bounded `[-1, 1]` window and choose
    // the least skewed valid one. This is a best-effort bounded search.
    for trans_corr in UNIMODULAR3_RANGE1.iter() {
        // trans_corr should keep centering translations
        if centering.lattice_points().iter().any(|translation| {
            let mut diff = trans_corr.linear.map(|e| e as f64) * translation - translation;
            diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
            diff.iter().any(|e| e.abs() > epsilon)
        }) {
            continue;
        }

        // trans_corr should keep the matrix representations of `conv_std_generators`.
        if conv_std_generators.iter().any(|ops| {
            let corr_ops = trans_corr.transform_operation(ops);
            if corr_ops.rotation != ops.rotation {
                true
            } else {
                let mut diff = corr_ops.translation - ops.translation;
                diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
                diff.iter().any(|e| e.abs() > epsilon)
            }
        }) {
            continue;
        }

        let refined_conv_trans = Transformation::new(
            transformation_to_prim_std.linear * centering.linear() * trans_corr.linear,
            transformation_to_prim_std.origin_shift,
        );
        let refined_conv_lattice = refined_conv_trans.transform_lattice(prim_lattice);
        // `skewness` gets smaller as the basis vectors are closer to orthorhombic.
        let skewness = refined_conv_lattice.lattice_constant()[3..]
            .iter()
            .map(|angle_deg| angle_deg.to_radians().cos().abs())
            .sum::<f64>();
        candidate_conv_transformations.push((skewness, centering.linear() * trans_corr.linear));
    }

    candidate_conv_transformations
        .into_iter()
        .min_by(|(skewness_lhs, _), (skewness_rhs, _)| {
            skewness_lhs.partial_cmp(skewness_rhs).unwrap()
        })
        .unwrap()
        .1
}

/// Test whether `position` matches a Wyckoff orbit described by `coordinates`
/// (the parsed-coordinates string from a Wyckoff entry, e.g. `"x,1/2,z"`)
/// modulo a bounded integer offset. Shared between the bulk and layer
/// pipelines: both Wyckoff databases store their representative coordinate as
/// the same `&str` shape, so the SNF-based offset search is identical.
///
/// The search probes offsets in `[-1, 1]^3` first, then the remaining shell
/// in `[-2, 2]^3`. The math (suppressing tolerances): we want `y, offset` so
/// that `lattice * (space.linear * y + space.origin - position - offset)` is
/// near zero. Decomposing `space.linear` as `D = L * space.linear * R` (SNF)
/// reduces the integer search to picking `D^{-1} L (offset + position -
/// space.origin)` and reading off whether the residual is small in Cartesian
/// distance.
pub(super) fn match_wyckoff_coordinates(
    position: &Position,
    coordinates: &str,
    lattice: &Lattice,
    symprec: f64,
) -> bool {
    let space = WyckoffPositionSpace::new(coordinates);
    let snf = SNF::new(&space.linear);

    let iter_multi_1 = iproduct!(-1..=1, -1..=1, -1..=1);
    let iter_multi_2 = iproduct!(-2_i32..=2_i32, -2_i32..=2_i32, -2_i32..=2_i32)
        .filter(|&(n1, n2, n3)| n1.abs() == 2 || n2.abs() == 2 || n3.abs() == 2);

    for offset in iter_multi_1.chain(iter_multi_2) {
        let offset = Vector3::new(offset.0 as f64, offset.1 as f64, offset.2 as f64);
        let b = snf.l.map(|e| e as f64) * (offset + position - space.origin);
        let mut rinvy = Vector3::zeros();
        for i in 0..3 {
            if snf.d[(i, i)] != 0 {
                rinvy[i] = b[i] / snf.d[(i, i)] as f64;
            }
        }

        let y = snf.r.map(|e| e as f64) * rinvy;
        let diff = space.linear.map(|e| e as f64) * y + space.origin - position - offset;
        if lattice.cartesian_coords(&diff).norm() < symprec {
            return true;
        }
    }
    false
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
