use super::point_group::{iter_trans_mat_basis, iter_unimodular_trans_mat};
use super::rotation_type::identify_rotation_type;
use super::space_group::match_origin_shift;
use crate::base::{project_rotations, Operations, UnimodularTransformation};

/// Return integral normalizer of the point group representative up to its centralizer.

/// Integral normalizers for all arithmetic point groups up to their centralizers.

/// Generate integral normalizer of the given space group up to its centralizer.
/// Because the factor group of the integral normalizer by the centralizer is isomorphic to a finite permutation group, the output is guaranteed to be finite.
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
