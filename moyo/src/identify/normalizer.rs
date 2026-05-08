use super::point_group::{iter_trans_mat_basis, iter_unimodular_trans_mat};
use super::rotation_type::identify_rotation_type;
use super::space_group::match_origin_shift;
use crate::base::{Operations, UnimodularTransformation, project_rotations};
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
