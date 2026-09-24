use nalgebra::Matrix3;
use std::collections::HashMap;

use super::conventional_coordinate_system::ConventionalCoordinateSystem;
use super::symmetrization::{symmetrize_lattice, symmetrize_positions};
use super::wyckoff::{assign_wyckoffs_by_orbit, group_sites_by_orbit, match_wyckoff_coordinates};
use crate::base::{
    Cell, MoyoError, Operations, Permutation, Transformation, UnimodularTransformation,
    project_rotations,
};
use crate::data::{HallNumber, WyckoffPosition, iter_wyckoff_positions};
use crate::identify::SpaceGroup;

/// Result of the full standardization pipeline for a primitive cell.
///
/// [`ConventionalCoordinateSystem`] selects the basis and origin. This type then
/// refines the lattice and atomic positions, constructs the primitive and conventional cells,
/// applies the requested Cartesian orientation, and assigns Wyckoff positions.
///
/// The [returned-cell specification](https://spglib.github.io/moyo/standardization/#returned-cell-specification)
/// describes the target symmetry, coordinate transformations, and Cartesian orientation.
///
/// Both output cells use the refined metric.
/// Coordinate transformations describe the selected basis and origin before
/// refinement; the rotation matrix describes only the final Cartesian rotation.
pub struct StandardizedCell {
    // ------------------------------------------------------------------------
    // Primitive standardized cell
    // ------------------------------------------------------------------------
    /// Primitive output, preserving the input primitive cell's site order and species.
    pub prim_cell: Cell,
    /// Coordinate transformation from the input primitive cell to the selected
    /// primitive system, before lattice and position refinement.
    pub prim_transformation: UnimodularTransformation,
    // ------------------------------------------------------------------------
    // Standardized cell
    // ------------------------------------------------------------------------
    /// Conventional output in the selected Hall setting.
    pub cell: Cell,
    /// One Wyckoff position per conventional site, in `cell` order and the selected setting.
    pub wyckoffs: Vec<WyckoffPosition>,
    /// Coordinate transformation from the input primitive cell to the selected
    /// conventional system, before lattice and position refinement.
    pub transformation: Transformation,
    // ------------------------------------------------------------------------
    // Miscellaneous
    // ------------------------------------------------------------------------
    /// Proper Cartesian rotation applied to both output lattices.
    /// Identity when `rotate_basis` is false; it does not encode lattice refinement.
    pub rotation_matrix: Matrix3<f64>,
    /// One primitive site index per conventional site: `site_mapping[j]` gives
    /// the site in `prim_cell` corresponding to site `j` in `cell`.
    pub site_mapping: Vec<usize>,
}

impl StandardizedCell {
    /// Standardize the input **primitive** cell.
    /// See [`Self`] for the returned-cell specification.
    ///
    /// Returns [`MoyoError::StandardizationError`] if lattice refinement fails.
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

        let centering = Transformation::from_linear(conv_trans_linear);
        let conv_lattice = centering.transform_lattice(&prim_std_cell_tmp.lattice);
        let (refined_lattice, rotation_matrix) =
            symmetrize_lattice(&conv_lattice, &project_rotations(&conv_std_operations))?;
        let (std_lattice, rotation_matrix) = if rotate_basis {
            (refined_lattice, rotation_matrix)
        } else {
            // Undo only the polar rotation, retaining the symmetric stretch.
            (
                refined_lattice.rotate(&rotation_matrix.transpose()),
                Matrix3::identity(),
            )
        };

        // Note: prim_transformation.transform_cell does not change the order of sites
        let prim_std_cell = Cell::new(
            centering.inverse_transform_lattice(&std_lattice),
            new_prim_std_positions,
            prim_std_cell_tmp.numbers.clone(),
        );

        // To (conventional) standardized cell
        let (std_cell, site_mapping) = centering.transform_cell(&prim_std_cell);

        Ok((
            prim_std_cell,
            prim_std_permutations,
            prim_transformation,
            std_cell,
            transformation,
            rotation_matrix,
            site_mapping,
        ))
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
