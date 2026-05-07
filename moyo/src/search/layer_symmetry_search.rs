use std::collections::{HashMap, HashSet, VecDeque};

use log::debug;

use super::layer_bravais_group::search_layer_bravais_group;
use super::solve::{
    PeriodicKdTree, pivot_site_indices, solve_correspondence,
    symmetrize_translation_from_permutation,
};
use crate::base::{
    AngleTolerance, Lattice, LayerCell, MoyoError, Operation, Operations, Permutation,
};

/// Coset representatives of the layer group of a primitive layer cell.
///
/// Mirrors [`super::PrimitiveSymmetrySearch`] but consumes the
/// LG-restricted Bravais group from [`search_layer_bravais_group`] (paper eq. 4
/// block form), so candidate rotations cannot mix the in-plane and aperiodic
/// directions.
#[derive(Debug)]
#[allow(dead_code)] // wired up by later layer-group milestones
pub struct LayerPrimitiveSymmetrySearch {
    /// Layer-group operations in the given primitive layer cell.
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
}

#[allow(dead_code)] // wired up by later layer-group milestones
impl LayerPrimitiveSymmetrySearch {
    /// Search layer-group coset representatives for a primitive layer cell.
    ///
    /// Assumes `primitive_layer_cell` is the output of
    /// [`super::LayerPrimitiveCell`] (its third basis vector equals the input
    /// `c` and its in-plane block is unreduced -- 2D Minkowski reduction is
    /// deferred to standardization).
    pub fn new(
        primitive_layer_cell: &LayerCell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        let inner_cell = primitive_layer_cell.cell();
        // Tolerance sanity check against the in-plane basis only -- `c` is not
        // a candidate translation direction for layer systems.
        let basis = &inner_cell.lattice.basis;
        let min_inplane_norm = basis.column(0).norm().min(basis.column(1).norm());
        let rough_symprec = 2.0 * symprec;
        if rough_symprec > min_inplane_norm / 2.0 {
            debug!(
                "symprec is too large compared to the in-plane basis vectors. Consider reducing symprec."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        let pkdtree = PeriodicKdTree::new(inner_cell, rough_symprec);
        let bravais_group =
            search_layer_bravais_group(&inner_cell.lattice, symprec, angle_tolerance)?;
        let pivot_site_indices = pivot_site_indices(&inner_cell.numbers);
        let mut symmetries_tmp = vec![];
        let src = pivot_site_indices[0];
        for rotation in bravais_group.iter() {
            let rotated_positions = inner_cell
                .positions
                .iter()
                .map(|pos| rotation.map(|e| e as f64) * pos)
                .collect::<Vec<_>>();
            for dst in pivot_site_indices.iter() {
                // Try to overlap the `src`-th site to the `dst`-th site.
                let translation = inner_cell.positions[*dst] - rotated_positions[src];
                let new_positions = rotated_positions
                    .iter()
                    .map(|pos| pos + translation)
                    .collect::<Vec<_>>();

                if let Some(permutation) =
                    solve_correspondence(&pkdtree, inner_cell, &new_positions)
                {
                    symmetries_tmp.push((*rotation, translation, permutation));
                    // WHY: keep all candidates -- multiple translations may match within rough tolerance.
                }
            }
        }
        debug!(
            "Number of layer-symmetry operation candidates: {}",
            symmetries_tmp.len()
        );

        // Tolerance on the fractional `c`-component of layer translations,
        // derived from the cartesian `symprec` divided by `|c|`.
        // WHY: layer-group operations have purely in-plane translations
        // (paper eq. 4 with t_3 = 0). A glide-mirror with non-zero t_z would
        // not be an element of any layer group.
        let nc = basis.column(2).norm();
        let z_tol = symprec / nc;

        // Purify symmetry operations by permutations and reject any with a
        // non-zero translation z-component.
        let mut operations_and_permutations = vec![];
        for (rotation, rough_translation, permutation) in symmetries_tmp.iter() {
            let (mut translation, distance) = symmetrize_translation_from_permutation(
                inner_cell,
                permutation,
                rotation,
                rough_translation,
            );
            if distance >= symprec {
                continue;
            }
            let mut tz = translation[2];
            tz -= tz.round();
            if tz.abs() > z_tol {
                continue;
            }
            translation[2] = 0.0;
            operations_and_permutations
                .push((Operation::new(*rotation, translation), permutation.clone()));
        }
        if operations_and_permutations.is_empty() {
            debug!(
                "No layer-symmetry operations are found. Consider increasing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooSmallToleranceError);
        }

        // Recover operations by group multiplication.
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        let mut operations = vec![];
        let mut permutations = vec![];
        queue.push_back((
            Operation::identity(),
            Permutation::identity(inner_cell.num_atoms()),
        ));
        while !queue.is_empty() {
            let (ops_lhs, permutation_lhs) = queue.pop_front().unwrap();
            if visited.contains(&ops_lhs.rotation) {
                continue;
            }
            visited.insert(ops_lhs.rotation);
            operations.push(ops_lhs.clone());
            permutations.push(permutation_lhs.clone());

            for (ops_rhs, permutation_rhs) in operations_and_permutations.iter() {
                let mut new_ops = ops_lhs.clone() * ops_rhs.clone();
                new_ops.translation -= new_ops.translation.map(|e| e.round());
                let new_permutation = permutation_lhs.clone() * permutation_rhs.clone();
                queue.push_back((new_ops, new_permutation));
            }
        }
        if operations.len() != operations_and_permutations.len() {
            debug!(
                "Found layer operations do not form a group. Consider reducing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        if !Self::check_closure(&operations, &inner_cell.lattice, rough_symprec) {
            debug!(
                "Some centering translations are missing. Consider reducing symprec and angle_tolerance."
            );
            return Err(MoyoError::TooLargeToleranceError);
        }

        debug!("Order of layer point group: {}", operations.len());
        Ok(Self {
            operations,
            permutations,
        })
    }

    fn check_closure(operations: &Operations, lattice: &Lattice, symprec: f64) -> bool {
        let mut translations_map = HashMap::new();
        for operation in operations.iter() {
            translations_map.insert(operation.rotation, operation.translation);
        }
        for ops1 in operations.iter() {
            for ops2 in operations.iter() {
                let ops12 = ops1.clone() * ops2.clone();
                let diff =
                    (translations_map[&ops12.rotation] - ops12.translation).map(|e| e - e.round());
                if lattice.cartesian_coords(&diff).norm() > symprec {
                    return false;
                }
            }
        }
        true
    }
}
