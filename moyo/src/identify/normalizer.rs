use itertools::Itertools;

use super::point_group::iter_trans_mat_basis;
use super::rotation_type::identify_rotation_type;
use super::space_group::match_origin_shift;
use crate::base::{project_rotations, Operations, UnimodularLinear, UnimodularTransformation};

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
        // Search integer linear combination such that the transformation matrix is unimodular
        // Consider coefficients in [-2, 2], which will be sufficient for Delaunay reduced basis
        for comb in (0..trans_mat_basis.len())
            .map(|_| -2..=2)
            .multi_cartesian_product()
        {
            let mut prim_trans_mat = UnimodularLinear::zeros();
            for (i, matrix) in trans_mat_basis.iter().enumerate() {
                prim_trans_mat += comb[i] * matrix;
            }
            let det = prim_trans_mat.map(|e| e as f64).determinant().round() as i32;
            if det < 0 {
                prim_trans_mat *= -1;
            }
            if det == 1 {
                if let Some(origin_shift) =
                    match_origin_shift(prim_operations, &prim_trans_mat, prim_generators, epsilon)
                {
                    conjugators.push(UnimodularTransformation::new(prim_trans_mat, origin_shift));
                }
                break;
            }
        }
    }
    conjugators
}
