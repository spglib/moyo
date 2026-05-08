use std::collections::HashMap;

use super::point_group::{iter_trans_mat_basis, iter_unimodular_trans_mat};
use super::rotation_type::identify_rotation_type;
use super::space_group::match_origin_shift;
use crate::base::{Operations, UnimodularLinear, UnimodularTransformation, project_rotations};
use crate::search::is_layer_block_form;

///
/// Return integral normalizer of the point group representative up to its centralizer.
///
/// Integral normalizers for all arithmetic point groups up to their centralizers.
/// Generate integral normalizer of the given space group up to its centralizer.
/// Because the factor group of the integral normalizer by the centralizer is isomorphic to a finite permutation group, the output is guaranteed to be finite.
///
/// Be careful that the input primitive operations should be in a reduced basis.
/// This function relies on the bounded search in `iter_unimodular_trans_mat`.
/// For a non-reduced basis, the output is only best effort and may not be exhaustive.
pub fn integral_normalizer(
    prim_operations: &Operations,
    prim_generators: &Operations,
    epsilon: f64,
) -> Vec<UnimodularTransformation> {
    let prim_rotations = project_rotations(prim_operations);
    let prim_rotation_generators = project_rotations(prim_generators);

    let rotation_types = prim_rotations
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();

    let mut conjugators = vec![];
    for trans_mat_basis in
        iter_trans_mat_basis(prim_rotations, rotation_types, prim_rotation_generators)
    {
        for prim_trans_mat in iter_unimodular_trans_mat(trans_mat_basis) {
            if let Some(origin_shift) =
                match_origin_shift(prim_operations, &prim_trans_mat, prim_generators, epsilon)
            {
                conjugators.push(UnimodularTransformation::new(prim_trans_mat, origin_shift));
                break;
            }
        }
    }
    conjugators
}

/// Layer-aware sibling of [`integral_normalizer`] for the (2 + 1)-dimensional
/// layer-group convention. Mirrors the bulk function with one extra filter:
/// only conjugators in **layer-block form** (`T[0,2] = T[1,2] = T[2,0] =
/// T[2,1] = 0`, `T[2,2] = +/-1`, paper Fu et al. 2024 eq. 4) are returned,
/// so the c-axis stays the aperiodic stacking direction.
///
/// Without the filter, `iter_unimodular_trans_mat`'s `[-2, 2]` enumeration
/// admits unimodular conjugators that rotate `c` into an in-plane axis
/// (e.g. for an oblique LG with a 2-fold along c, any 3D `(P, p)`
/// permuting c with a or b is in the bulk integral normalizer of the
/// rotation set). Such conjugators are valid as bulk normalizer elements
/// but break the periodicity contract enforced by `LayerLattice` /
/// `LayerCell`. The `_2_1` suffix denotes the (2 + 1)-D convention; layer
/// cells are still 3D structures, so `_2d` would be misleading.
///
/// Same input contract as [`integral_normalizer`]: `prim_operations` is
/// expected in a reduced basis; the result is best-effort within
/// `iter_unimodular_trans_mat`'s bounded search.
pub fn integral_normalizer_2_1(
    prim_operations: &Operations,
    prim_generators: &Operations,
    epsilon: f64,
) -> Vec<UnimodularTransformation> {
    let prim_rotations = project_rotations(prim_operations);
    let prim_rotation_generators = project_rotations(prim_generators);

    let rotation_types = prim_rotations
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();

    let mut conjugators = vec![];
    for trans_mat_basis in
        iter_trans_mat_basis(prim_rotations, rotation_types, prim_rotation_generators)
    {
        for prim_trans_mat in iter_unimodular_trans_mat(trans_mat_basis) {
            if !is_layer_block_form(&prim_trans_mat) {
                continue;
            }
            if let Some(mut origin_shift) =
                match_origin_shift(prim_operations, &prim_trans_mat, prim_generators, epsilon)
            {
                // `match_origin_shift` -> `solve_mod1` reduces mod 1 on all
                // three components. The layer cell's c-axis is aperiodic, so
                // override the c-component with the exact value derived from
                // a c-flipping generator's c-row -- see `layer_exact_s_z`
                // for the derivation, and `docs/layer_architecture.md` for
                // why the override is baked in here rather than at the
                // caller.
                if let Some(exact_s_z) =
                    layer_exact_s_z(prim_operations, &prim_trans_mat, prim_generators, epsilon)
                {
                    origin_shift[2] = prim_trans_mat[(2, 2)] as f64 * exact_s_z;
                }
                conjugators.push(UnimodularTransformation::new(prim_trans_mat, origin_shift));
                break;
            }
        }
    }
    conjugators
}

/// Recover the exact `s_z` (origin-shift c-component) by reading any
/// c-flipping generator's c-row of the matching equation.
///
/// Returns `None` when no c-flipping generator exists in the input ops or
/// the database (e.g. LG 1-5 oblique without inversion / horizontal
/// mirror), in which case `s_z` is genuinely unconstrained and the modular
/// solver's `s_z = 0` answer is correct.
///
/// When multiple c-flipping generators exist a debug-assert checks that
/// they yield the same `s_z` (a group-consistency invariant).
fn layer_exact_s_z(
    prim_layer_operations: &Operations,
    trans_mat: &UnimodularLinear,
    db_prim_generators: &Operations,
    epsilon: f64,
) -> Option<f64> {
    // Mirror the front of `match_origin_shift`: apply `trans_mat` to the
    // input ops and index by rotation so we can look up `t_target` per
    // generator.
    let new_prim_operations = UnimodularTransformation::from_linear(*trans_mat)
        .transform_operations(prim_layer_operations);
    let mut hm_translations = HashMap::new();
    for op in new_prim_operations.iter() {
        hm_translations.insert(op.rotation, op.translation);
    }

    let mut chosen: Option<f64> = None;
    for db_op in db_prim_generators.iter() {
        if db_op.rotation[(2, 2)] != -1 {
            continue;
        }
        let target_t = hm_translations.get(&db_op.rotation)?;
        // (W - E) c-row for any layer-block W with W[2,2] = -1 is
        // (0, 0, -2); so -2 s_z = t_db_z - t_target_z, i.e.
        // s_z = (t_target_z - t_db_z) / 2.
        let s_z = (target_t[2] - db_op.translation[2]) / 2.0;
        match chosen {
            None => chosen = Some(s_z),
            Some(prev) => debug_assert!(
                (prev - s_z - (prev - s_z).round()).abs() < epsilon.max(1e-6),
                "c-flipping generators disagree on s_z (prev={}, new={})",
                prev,
                s_z
            ),
        }
    }
    chosen
}
