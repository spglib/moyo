use super::layer_bravais_group::is_layer_block_form;
use super::primitive_symmetry_search::PrimitiveSymmetrySearch;
use crate::base::{
    AngleTolerance, Cell, Lattice, LayerCell, MoyoError, Operation, Operations, Permutation,
};

/// Coset representatives of the layer group of a primitive layer cell.
///
/// Wraps `PrimitiveSymmetrySearch::new` and post-filters its operations to
/// the layer-group ones: rotations must satisfy paper Fu et al. 2024 eq. 4
/// (block form, in/out separated; see `is_layer_block_form`) and the
/// intrinsic c-translation must vanish. The intrinsic c-translation is the
/// component of `t` along the c invariant subspace of `W`: for `W[2,2] = +1`
/// it is `t_3` itself (so a non-zero `t_3` is a c-axis screw, forbidden in
/// a layer group), while for `W[2,2] = -1` it averages to zero by
/// construction (`(t_3 + (-t_3))/2 = 0`), leaving `t_3` as a free origin
/// shift along the aperiodic axis. The resulting subset is a subgroup of
/// the bulk space group, hence still closed; no separate closure check is
/// needed here.
#[derive(Debug)]
pub struct LayerPrimitiveSymmetrySearch {
    /// Layer-group operations in the given primitive layer cell.
    pub operations: Operations,
    pub permutations: Vec<Permutation>,
}

impl LayerPrimitiveSymmetrySearch {
    /// Search layer-group coset representatives for a primitive layer cell.
    ///
    /// Assumes `primitive_layer_cell` is the output of
    /// `super::LayerPrimitiveCell` (its third basis vector equals the input
    /// `c` and its in-plane block is unreduced -- 2D Minkowski reduction is
    /// deferred to standardization).
    pub fn new(
        primitive_layer_cell: &LayerCell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
    ) -> Result<Self, MoyoError> {
        // Reconstruct a bulk `Cell` once: explicit at the call site so the
        // layer-to-bulk crossing is visible.
        let cell = Cell::new(
            Lattice {
                basis: *primitive_layer_cell.lattice().basis(),
            },
            primitive_layer_cell.positions().to_vec(),
            primitive_layer_cell.numbers().to_vec(),
        );

        // Reuse the bulk space-group search. For a layer-shaped input the
        // candidate bulk Bravais group is a superset of the layer one; the
        // post-filter below prunes the non-layer elements.
        let bulk = PrimitiveSymmetrySearch::new(&cell, symprec, angle_tolerance)?;

        let nc = primitive_layer_cell.lattice().basis().column(2).norm();
        let z_tol = symprec / nc;

        let mut operations = vec![];
        let mut permutations = vec![];
        for (op, perm) in bulk.operations.iter().zip(bulk.permutations.iter()) {
            if !is_layer_block_form(&op.rotation) {
                continue;
            }
            // Reject c-axis screws: when W[2,2] = +1 the intrinsic z-translation
            // is t_3 itself and must vanish. For W[2,2] = -1 the intrinsic part
            // is automatically zero, so t_3 is just the location of the
            // symmetry center along the aperiodic c direction.
            if op.rotation[(2, 2)] == 1 {
                let mut tz = op.translation[2];
                tz -= tz.round();
                if tz.abs() > z_tol {
                    continue;
                }
            }
            operations.push(Operation::new(op.rotation, op.translation));
            permutations.push(perm.clone());
        }

        if operations.is_empty() {
            return Err(MoyoError::TooSmallToleranceError);
        }
        Ok(Self {
            operations,
            permutations,
        })
    }
}
