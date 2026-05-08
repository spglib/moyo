use super::layer_bravais_group::is_layer_block_form;
use super::primitive_symmetry_search::PrimitiveSymmetrySearch;
use crate::base::{
    AngleTolerance, Cell, Lattice, LayerCell, MoyoError, Operation, Operations, Permutation,
    Rotation,
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
            // Recover the exact c-component of the translation from atom
            // images: the bulk search returns `t` reduced mod 1 in all three
            // axes, but for layer cells `c` is aperiodic and the inversion-
            // center / mirror-plane location encoded in `t_z` is lost when
            // mod-1 reduction collapses centers at z and z + 1/2 into the
            // same translation. For layer-block W (`W[2,*] = (0, 0, +/-1)`),
            // `(W r)_z = W[2,2] * r_z`, so the exact t_z is the orbit-mean
            // `r_image_z - W[2,2] * r_z`.
            let exact_tz = exact_tz(primitive_layer_cell, &op.rotation, perm);
            // Reject c-axis screws: when W[2,2] = +1 the intrinsic z-translation
            // is t_3 itself and must vanish. For W[2,2] = -1 the intrinsic part
            // is automatically zero, so t_3 is just the location of the
            // symmetry center along the aperiodic c direction. The screw check
            // operates on `exact_tz` modulo lattice (an exact `t_z = 1.0` is a
            // unit-cell translation, not a screw), so the `tz - tz.round()`
            // reduction is still correct after the lift.
            if op.rotation[(2, 2)] == 1 {
                let tz = exact_tz - exact_tz.round();
                if tz.abs() > z_tol {
                    continue;
                }
            }
            let mut translation = op.translation;
            translation[2] = exact_tz;
            operations.push(Operation::new(op.rotation, translation));
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

/// Recover the exact c-component of an operation's translation in a
/// `LayerCell`. For layer-block `W` (`W[2,*] = (0, 0, W[2,2])`,
/// `W[2,2] = +/- 1`), `(W r)_z = W[2,2] * r_z`, so the displacement
/// `r_image_z - W[2,2] * r_z` is the orbit-invariant exact c-translation.
/// Unlike the in-plane components, no `mod 1` reduction is applied: `c`
/// is aperiodic and the raw value carries the inversion-center /
/// mirror-plane location.
///
/// Atom positions are stored in `[0, 1)` even for the aperiodic c-axis,
/// so an atom at "real" `z = -0.04` is stored at `z = 0.96`. This wrap
/// inflates the per-atom contribution by an integer for atoms whose
/// image crosses a cell boundary. We anchor to atom 0 and round each
/// subsequent atom's contribution to atom 0's reference, then average.
fn exact_tz(cell: &LayerCell, rotation: &Rotation, permutation: &Permutation) -> f64 {
    let n = cell.num_atoms();
    debug_assert!(n > 0);
    let w22 = rotation[(2, 2)] as f64;
    let positions = cell.positions();
    let raw_tz = |i: usize| positions[permutation.apply(i)][2] - w22 * positions[i][2];
    let reference = raw_tz(0);
    let mut sum = reference;
    for i in 1..n {
        let candidate = raw_tz(i);
        sum += candidate - (candidate - reference).round();
    }
    sum / n as f64
}
