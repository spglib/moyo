use itertools::Itertools;
use once_cell::sync::Lazy;
use std::collections::HashSet;

use super::rotation_type::identify_rotation_type;
use crate::base::{traverse, Rotations, UnimodularLinear};
use crate::data::{ArithmeticNumber, PointGroupRepresentative};
use crate::math::sylvester3;

/// Return integral normalizer of the point group representative up to its centralizer.
pub fn get_point_group_normalizer(
    arithmetic_number: ArithmeticNumber,
) -> Option<Vec<UnimodularLinear>> {
    POINT_GROUP_NORMALIZERS
        .get((arithmetic_number - 1) as usize)
        .cloned()
}

/// Integral normalizers for all arithmetic point groups up to their centralizers.
static POINT_GROUP_NORMALIZERS: Lazy<Vec<Vec<UnimodularLinear>>> = Lazy::new(|| {
    (1..=73)
        .map(|arithmetic_number| {
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
            let prim_generators = point_group_db.primitive_generators();
            let prim_rotations = traverse(&prim_generators);

            let mut prim_rotations_set = HashSet::new();
            for rotation in prim_rotations.iter() {
                prim_rotations_set.insert(rotation.clone());
            }
            integral_normalizer(&prim_rotations, &prim_generators)
        })
        .collect()
});

#[allow(dead_code)]
/// Generate integral normalizer of the given point group up to its centralizer.
/// Because the factor group of the integral normalizer by the centralizer is isomorphic to a finite permutation group, the output is guaranteed to be finite.
fn integral_normalizer(
    prim_rotations: &Rotations,
    prim_generators: &Rotations,
) -> Vec<UnimodularLinear> {
    let rotation_types = prim_rotations
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();
    let prim_generator_rotation_types = prim_generators
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();

    // Try to map generators
    let order = prim_rotations.len();
    let candidates: Vec<Vec<usize>> = prim_generator_rotation_types
        .iter()
        .map(|&rotation_type| {
            (0..order)
                .filter(|&i| rotation_types[i] == rotation_type)
                .collect()
        })
        .collect();

    // TODO: unify with match_with_point_group
    let mut conjugators = vec![];
    for pivot in candidates
        .iter()
        .map(|v| v.iter())
        .multi_cartesian_product()
    {
        // Solve P^-1 * prim_rotations[pivot[i]] * P = prim_generators[i] (for all i)
        if let Some(trans_mat_basis) = sylvester3(
            &pivot
                .iter()
                .map(|&i| prim_rotations[*i])
                .collect::<Vec<_>>(),
            prim_generators,
        ) {
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
                    conjugators.push(prim_trans_mat);
                    break;
                }
            }
        }
    }
    conjugators
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::base::traverse;
    use crate::data::PointGroupRepresentative;

    #[test]
    fn test_integral_normalizer() {
        for arithmetic_number in 1..=73 {
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
            let prim_generators = point_group_db.primitive_generators();
            let prim_rotations = traverse(&prim_generators);

            let mut prim_rotations_set = HashSet::new();
            for rotation in prim_rotations.iter() {
                prim_rotations_set.insert(rotation.clone());
            }

            let normalizer = get_point_group_normalizer(arithmetic_number).unwrap();
            assert!(normalizer.len() > 0);

            for linear in normalizer.iter() {
                let linear_inv = linear
                    .map(|e| e as f64)
                    .try_inverse()
                    .unwrap()
                    .map(|e| e.round() as i32);
                let prim_rotations_actual: Vec<_> = prim_rotations
                    .iter()
                    .map(|r| linear_inv * r * linear)
                    .collect();
                for rotation in prim_rotations_actual {
                    assert!(prim_rotations_set.contains(&rotation));
                }
            }
        }
    }
}
