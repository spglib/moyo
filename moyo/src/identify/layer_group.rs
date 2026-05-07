use log::debug;
use serde::Serialize;

use super::point_group::PointGroup;
use super::space_group::match_origin_shift;
use crate::base::{
    MoyoError, Operations, UnimodularLinear, UnimodularTransformation, project_rotations,
};
use crate::data::{
    LayerHallNumber, LayerHallSymbol, LayerNumber, LayerSetting, arithmetic_crystal_class_entry,
    layer_arithmetic_crystal_class_entry, layer_hall_symbol_entry,
};

/// Identified layer-group type for a primitive layer cell.
///
/// Mirrors [`super::SpaceGroup`] for the bulk space-group case: result
/// carries the layer-group number (1..=80, paper Fu et al. 2024 Table 4)
/// and a representative [`LayerHallNumber`] (1..=116) along with the
/// `UnimodularTransformation` that maps the input primitive layer cell to
/// the database canonical for that Hall.
///
/// [`LayerGroup::new`] runs the bulk [`PointGroup::new`] to extract the
/// geometric crystal class, uses it to pre-filter Hall candidates, and for
/// each survivor sweeps a small set of layer-valid normalizer corrections
/// (`identity` and an in-plane `a` ↔ `b` swap) before attempting
/// `match_origin_shift`. This handles the common axis-labelling
/// ambiguities -- monoclinic-rectangular `:a` vs `:b` and orthorhombic
/// `abc` vs `b-ac` -- when the input lies in a non-canonical basis.
#[derive(Debug, Clone, Serialize)]
pub struct LayerGroup {
    pub number: LayerNumber,
    pub hall_number: LayerHallNumber,
    /// Transformation from the input primitive basis to the canonical
    /// primitive basis of `hall_number`.
    pub transformation: UnimodularTransformation,
}

impl LayerGroup {
    /// Identify a layer group from primitive layer-cell operations.
    pub fn new(
        prim_layer_operations: &Operations,
        setting: LayerSetting,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        // Identify the bulk point group first (mirrors `SpaceGroup::new`).
        // Only the geometric crystal class is consumed: the bulk arithmetic
        // class would normalise to a 3D-canonical orientation that does
        // not preserve the layer's c-axis identity, so its `prim_trans_mat`
        // is not fed into `match_origin_shift` (every layer Hall in the
        // database is stored in the layer-canonical basis with c aperiodic).
        let prim_rotations = project_rotations(prim_layer_operations);
        let point_group = PointGroup::new(&prim_rotations)?;
        let geometric_crystal_class = arithmetic_crystal_class_entry(point_group.arithmetic_number)
            .ok_or(MoyoError::LayerGroupTypeIdentificationError)?
            .geometric_crystal_class;
        debug!(
            "Layer point group: geometric crystal class {:?}",
            geometric_crystal_class
        );

        let corrections = layer_correction_matrices();

        for hall_number in setting.hall_numbers() {
            let entry = layer_hall_symbol_entry(hall_number)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
            let layer_arith = layer_arithmetic_crystal_class_entry(entry.arithmetic_number)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
            if layer_arith.geometric_crystal_class != geometric_crystal_class {
                continue;
            }

            let lh_symbol = LayerHallSymbol::from_hall_number(hall_number)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
            let db_prim_generators = lh_symbol.primitive_generators();

            for trans_mat in &corrections {
                if let Some(origin_shift) = match_origin_shift(
                    prim_layer_operations,
                    trans_mat,
                    &db_prim_generators,
                    epsilon,
                ) {
                    debug!(
                        "Matched layer Hall number {} (LG {}) via correction {:?}",
                        hall_number, entry.number, trans_mat
                    );
                    return Ok(Self {
                        number: entry.number,
                        hall_number,
                        transformation: UnimodularTransformation::new(*trans_mat, origin_shift),
                    });
                }
            }
        }
        Err(MoyoError::LayerGroupTypeIdentificationError)
    }
}

/// Layer-allowed normalizer transformations for the basis search in
/// [`LayerGroup::new`]. Both elements preserve the layer block form
/// `W_i3 = W_3i = 0, W_33 = ±1` (paper Fu et al. 2024 eq. 4), so feeding
/// them into [`match_origin_shift`] does not break the c-aperiodic
/// invariant.
///
/// - `identity`: input already in the canonical layer basis.
/// - `b a -c`: in-plane `a` ↔ `b` swap, paired with a c sign-flip to keep
///   the determinant equal to +1. Covers the monoclinic-rectangular
///   `:a` vs `:b` axis labelling and the orthorhombic `abc` vs `b-ac`
///   axis swap. For tetragonal / hexagonal point groups the swap is
///   already an element of the point group, so the second candidate is
///   redundant but harmless (it will fall through to the next Hall).
fn layer_correction_matrices() -> [UnimodularLinear; 2] {
    [
        UnimodularLinear::identity(),
        UnimodularLinear::new(
            0, 1, 0, //
            1, 0, 0, //
            0, 0, -1, //
        ),
    ]
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use nalgebra::vector;

    use super::*;
    use crate::data::{LayerHallSymbol, iter_layer_hall_symbol_entry};

    /// Round-trip: for every Hall number in LAYER_HALL_SYMBOL_DATABASE,
    /// expand its symbol to primitive operations and pass them through
    /// `LayerGroup::new`; the returned `(number, hall_number)` must
    /// match the database entry, and the identified transformation must
    /// be unimodular and recover the same operations modulo translation.
    #[test]
    fn test_round_trip_all_layer_hall_numbers() {
        for entry in iter_layer_hall_symbol_entry() {
            let lh_symbol = LayerHallSymbol::from_hall_number(entry.hall_number).unwrap();
            let prim_operations = lh_symbol.primitive_traverse();

            let identified = LayerGroup::new(
                &prim_operations,
                LayerSetting::HallNumber(entry.hall_number),
                1e-8,
            )
            .unwrap_or_else(|_| {
                panic!(
                    "failed to identify layer hall_number {} (LG {})",
                    entry.hall_number, entry.number
                )
            });

            assert_eq!(
                identified.hall_number, entry.hall_number,
                "Hall number mismatch for LG {}",
                entry.number
            );
            assert_eq!(
                identified.number, entry.number,
                "LG number mismatch for hall_number {}",
                entry.hall_number
            );
            assert_eq!(
                identified
                    .transformation
                    .linear_as_f64()
                    .determinant()
                    .round() as i32,
                1,
                "non-unimodular transformation for hall_number {}",
                entry.hall_number
            );

            // Verify the transformation maps input operations onto the
            // database operations modulo unit-cell translations.
            let matched_lh = LayerHallSymbol::from_hall_number(identified.hall_number).unwrap();
            let matched_prim_operations = matched_lh.primitive_traverse();
            let mut hm_translations = HashMap::new();
            for op in matched_prim_operations.iter() {
                hm_translations.insert(op.rotation, op.translation);
            }

            let transformed = identified
                .transformation
                .transform_operations(&prim_operations);
            assert_eq!(
                matched_prim_operations.len(),
                transformed.len(),
                "operation count mismatch for hall_number {}",
                entry.hall_number
            );
            for op in transformed.iter() {
                let target = hm_translations.get(&op.rotation).unwrap_or_else(|| {
                    panic!(
                        "rotation missing in DB for hall_number {}",
                        entry.hall_number
                    )
                });
                let mut diff = *target - op.translation;
                diff -= diff.map(|e| e.round());
                assert_relative_eq!(diff, vector![0.0, 0.0, 0.0], epsilon = 1e-8);
            }
        }
    }

    /// With `LayerSetting::Standard`, an input expanded from the LG's
    /// default Hall must round-trip to the same Hall number (sanity check
    /// that the default-iteration path works).
    #[test]
    fn test_round_trip_default_hall_numbers_via_standard_setting() {
        let standard_halls = LayerSetting::Standard.hall_numbers();
        for &hall_number in &standard_halls {
            let lh_symbol = LayerHallSymbol::from_hall_number(hall_number).unwrap();
            let prim_operations = lh_symbol.primitive_traverse();
            let identified =
                LayerGroup::new(&prim_operations, LayerSetting::Standard, 1e-8).unwrap();
            assert_eq!(identified.hall_number, hall_number);
        }
    }

    /// Inputs in a non-canonical basis (the `:b` axis labelling for
    /// monoclinic-rectangular, the `b-ac` axis swap for orthorhombic)
    /// should round-trip to the LG's Standard Hall via the in-plane
    /// `a` ↔ `b` correction matrix in [`layer_correction_matrices`].
    /// The reported Hall number is the Standard one (i.e. the canonical
    /// `:a` / `abc` row), not the input Hall, since identification chose
    /// the canonical-basis match.
    ///
    /// Restricted to **`p`-centered** Halls. For `c`-centered LGs the
    /// conjugation of `b a -c` by the centering transform reduces to a
    /// group-internal rotation in the primitive basis, so the
    /// non-canonical-basis alignment for those LGs needs a different
    /// (centering-aware) correction matrix; tracked as future work.
    #[test]
    fn test_normalizer_aligns_non_canonical_basis() {
        let cases = [
            // (input Hall, expected LG number, expected Standard Hall)
            // Monoclinic-rectangular :b -> :a (p-centered).
            (13, 8, 12),  // LG 8:  p 1 2 1 -> p 2 1 1
            (15, 9, 14),  // LG 9:  p 1 2_1 1 -> p 2_1 1 1
            (19, 11, 18), // LG 11: p 1 m 1 -> p m 1 1
            (21, 12, 20), // LG 12: p 1 a 1 -> p b 1 1
            (25, 14, 24), // LG 14: p 1 2/m 1 -> p 2/m 1 1
            (27, 15, 26), // LG 15: p 1 2_1/m 1 -> p 2_1/m 1 1
            (29, 16, 28), // LG 16: p 1 2/a 1 -> p 2/b 1 1
            (31, 17, 30), // LG 17: p 1 2_1/a 1 -> p 2_1/b 1 1
            // Orthorhombic abc <-> b-ac (p-centered).
            (36, 20, 35), // LG 20: p 2 2_1 2 -> p 2_1 2 2
            (41, 24, 40), // LG 24: p b m 2 -> p m a 2
            (45, 27, 44), // LG 27: p 2 m m -> p m 2 m
            (66, 38, 65), // LG 38: p b m b -> p m a a
        ];
        for (input_hall, expected_lg, expected_standard_hall) in cases {
            let lh_symbol = LayerHallSymbol::from_hall_number(input_hall).unwrap();
            let prim_operations = lh_symbol.primitive_traverse();
            let identified = LayerGroup::new(&prim_operations, LayerSetting::Standard, 1e-8)
                .unwrap_or_else(|_| panic!("non-canonical-basis case hall={} failed", input_hall));
            assert_eq!(
                identified.number, expected_lg,
                "input hall {} (LG {}): wrong LG number",
                input_hall, expected_lg
            );
            assert_eq!(
                identified.hall_number, expected_standard_hall,
                "input hall {} (LG {}): expected to land on Standard hall {} after correction",
                input_hall, expected_lg, expected_standard_hall
            );
            // The recovered transformation must be unimodular.
            assert_eq!(
                identified
                    .transformation
                    .linear_as_f64()
                    .determinant()
                    .round() as i32,
                1,
                "non-unimodular transformation for input hall {}",
                input_hall
            );
        }
    }
}
