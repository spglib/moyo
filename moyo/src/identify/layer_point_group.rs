use log::debug;
use serde::{Deserialize, Serialize};

use super::point_group::{
    identify_geometric_crystal_class, iter_trans_mat_basis, iter_unimodular_trans_mat,
};
use super::rotation_type::identify_rotation_type;
use crate::base::{MoyoError, Rotations, UnimodularLinear};
use crate::data::{
    LayerArithmeticNumber, LayerPointGroupRepresentative, iter_layer_arithmetic_crystal_entry,
};
use crate::search::is_layer_block_form;

/// Crystallographic layer point group with arithmetic-class identification.
///
/// Mirrors [`super::PointGroup`] for the bulk space-group case but iterates
/// the **layer** arithmetic classes (1..=43, paper Fu et al. 2024 Table 4)
/// with canonical generators in the layer convention (the aperiodic axis
/// is always `c`).
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct LayerPointGroup {
    pub arithmetic_number: LayerArithmeticNumber,
    /// Layer-block-form unimodular matrix taking the input primitive
    /// basis to the canonical primitive basis of `arithmetic_number`.
    pub prim_trans_mat: UnimodularLinear,
}

impl LayerPointGroup {
    /// Identify the layer arithmetic crystal class from primitive
    /// layer-block rotations.
    ///
    /// `prim_rotations` are assumed to be in layer-block form already
    /// (the layer search step filters bulk operations through
    /// `is_layer_block_form` before exposing them); this routine does
    /// not re-validate.
    pub fn new(prim_rotations: &Rotations) -> Result<Self, MoyoError> {
        let rotation_types = prim_rotations.iter().map(identify_rotation_type).collect();
        let geometric_crystal_class = identify_geometric_crystal_class(&rotation_types)?;
        debug!(
            "Layer point group: geometric crystal class {:?}",
            geometric_crystal_class
        );

        for entry in iter_layer_arithmetic_crystal_entry() {
            if entry.geometric_crystal_class != geometric_crystal_class {
                continue;
            }

            let point_group_db = LayerPointGroupRepresentative::from_layer_arithmetic_crystal_class(
                entry.arithmetic_number,
            );
            let prim_generators_db = point_group_db.primitive_generators();

            // Try identity first.
            if prim_generators_db
                .iter()
                .all(|r| prim_rotations.contains(r))
            {
                return Ok(LayerPointGroup {
                    arithmetic_number: entry.arithmetic_number,
                    prim_trans_mat: UnimodularLinear::identity(),
                });
            }

            // Sylvester search for a layer-block-form conjugator.
            for trans_mat_basis in iter_trans_mat_basis(
                prim_rotations.clone(),
                rotation_types.clone(),
                prim_generators_db,
            ) {
                for prim_trans_mat in iter_unimodular_trans_mat(trans_mat_basis) {
                    // Filter for layer-block form: an unconstrained
                    // linear combination of mixed-form Sylvester basis
                    // matrices may produce a unimodular `T` that
                    // rotates `c` into an in-plane axis. Such `T` is
                    // valid as a 3D unimodular transformation but
                    // breaks the periodicity contract enforced by
                    // `LayerLattice` / `LayerCell`.
                    if !is_layer_block_form(&prim_trans_mat) {
                        continue;
                    }
                    return Ok(LayerPointGroup {
                        arithmetic_number: entry.arithmetic_number,
                        prim_trans_mat,
                    });
                }
            }
        }

        Err(MoyoError::ArithmeticCrystalClassIdentificationError)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::base::traverse;

    /// Round-trip: for every layer arithmetic class, expand its
    /// canonical primitive generators, traverse to the full rotation
    /// group, and confirm `LayerPointGroup::new` returns the same
    /// arithmetic_number with a unimodular layer-block prim_trans_mat.
    #[test]
    fn test_layer_point_group_round_trip_all_arithmetic_classes() {
        for arithmetic_number in 1..=43 {
            let rep = LayerPointGroupRepresentative::from_layer_arithmetic_crystal_class(
                arithmetic_number,
            );
            let prim_generators = rep.primitive_generators();
            let prim_rotations = traverse(&prim_generators);

            let identified = LayerPointGroup::new(&prim_rotations).unwrap_or_else(|_| {
                panic!(
                    "LayerPointGroup::new failed on canonical input for arithmetic_number={}",
                    arithmetic_number
                )
            });

            assert_eq!(
                identified.arithmetic_number, arithmetic_number,
                "arithmetic_number mismatch for canonical input arithmetic_number={}",
                arithmetic_number
            );
            assert_eq!(
                identified
                    .prim_trans_mat
                    .map(|e| e as f64)
                    .determinant()
                    .round() as i32,
                1,
                "non-unimodular prim_trans_mat for arithmetic_number={}",
                arithmetic_number
            );
            assert!(
                is_layer_block_form(&identified.prim_trans_mat),
                "non-layer-block prim_trans_mat for arithmetic_number={}: {:?}",
                arithmetic_number,
                identified.prim_trans_mat
            );

            // Conjugate the input rotations by prim_trans_mat and check
            // they cover the canonical primitive generators (modulo
            // group element relabeling).
            let prim_trans_inv = identified
                .prim_trans_mat
                .map(|e| e as f64)
                .try_inverse()
                .unwrap()
                .map(|e| e.round() as i32);
            let conjugated: HashSet<_> = prim_rotations
                .iter()
                .map(|r| prim_trans_inv * r * identified.prim_trans_mat)
                .collect();
            for g in &prim_generators {
                assert!(
                    conjugated.contains(g),
                    "conjugated rotation set does not cover db generator {:?} \
                     for arithmetic_number={}",
                    g,
                    arithmetic_number
                );
            }
        }
    }
}
