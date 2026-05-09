use log::debug;
use serde::Serialize;

use super::layer_point_group::LayerPointGroup;
use super::normalizer::integral_normalizer_2_1;
use crate::base::{MoyoError, Operations, UnimodularTransformation, project_rotations};
use crate::data::{
    LayerHallNumber, LayerHallSymbol, LayerNumber, LayerSetting, layer_hall_symbol_entry,
};

/// Identified layer-group type for a primitive layer cell.
///
/// Mirrors [`crate::identify::SpaceGroup`] for the bulk space-group case: result
/// carries the layer-group number (1..=80, paper Fu et al. 2024 Table 4)
/// and a representative [`LayerHallNumber`] (1..=116) along with the
/// `UnimodularTransformation` that maps the input primitive layer cell to
/// the database canonical for that Hall.
///
/// [`LayerGroup::new`] is the two-stage pipeline described in
/// `docs/layer_architecture.md`:
///
/// 1. [`super::LayerPointGroup::new`] matches the layer arithmetic
///    crystal class (1..=43) and synthesises a layer-block-form
///    `prim_trans_mat` that aligns the input primitive basis with the
///    canonical primitive basis of the matched class.
/// 2. [`super::integral_normalizer_2_1`] enumerates within-class
///    layer-block conjugators (axis settings, monoclinic cell choices,
///    C-centered shears, the trigonal axis-flip) using the same
///    Sylvester + `match_origin_shift` flow that the bulk
///    `MagneticSpaceGroup::new` Type-III branch uses, and bakes the
///    aperiodic-c origin-shift correction into each returned
///    conjugator. For each candidate Hall in `setting.hall_numbers()`
///    whose arithmetic class matches Stage 1, the first conjugator
///    returned is the answer.
///
/// The static `LAYER_CORRECTION_MATRICES` constant from earlier
/// iterations is gone: every correction it covered (identity,
/// `b a -c`, `diag(-1, +1, -1)`) is a specific element of the
/// integral normalizer Stage 2 enumerates, plus Stage 2 covers the
/// C-centered shear cases that no finite static list could.
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
        // Stage 1: match the layer arithmetic crystal class.
        let prim_rotations = project_rotations(prim_layer_operations);
        let layer_point_group = LayerPointGroup::new(&prim_rotations)?;
        debug!(
            "Layer point group: arithmetic_number={}",
            layer_point_group.arithmetic_number
        );

        // Stages 2 + 3: for each candidate Hall, enumerate layer-block
        // normalizer conjugators that complete the SG-level match, and
        // post-fix the aperiodic-c branch ambiguity.
        for hall_number in setting.hall_numbers() {
            let entry = layer_hall_symbol_entry(hall_number).unwrap();
            if entry.arithmetic_number != layer_point_group.arithmetic_number {
                continue;
            }

            let lh_symbol = LayerHallSymbol::from_hall_number(hall_number)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
            let db_prim_generators = lh_symbol.primitive_generators();

            let conjugators =
                integral_normalizer_2_1(prim_layer_operations, &db_prim_generators, epsilon);
            let Some(transformation) = conjugators.into_iter().next() else {
                continue;
            };
            debug!(
                "Matched layer Hall number {} (LG {}) via prim_trans_mat {:?}",
                hall_number, entry.number, transformation.linear
            );
            return Ok(Self {
                number: entry.number,
                hall_number,
                transformation,
            });
        }
        Err(MoyoError::LayerGroupTypeIdentificationError)
    }
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
    /// should round-trip to the LG's Standard Hall via the integral
    /// normalizer enumeration. The reported Hall number is the Standard
    /// one (i.e. the canonical `:a` / `abc` row), not the input Hall,
    /// since identification chose the canonical-basis match.
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
