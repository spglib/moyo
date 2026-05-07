use log::debug;
use serde::Serialize;

use super::space_group::match_origin_shift;
use crate::base::{MoyoError, Operations, UnimodularLinear, UnimodularTransformation};
use crate::data::{
    LayerHallNumber, LayerHallSymbol, LayerNumber, LayerSetting, layer_hall_symbol_entry,
};

/// Identified layer-group type for a primitive layer cell.
///
/// Mirrors [`super::SpaceGroup`] for the bulk space-group case: result
/// carries the layer-group number (1..=80, paper Fu et al. 2024 Table 4)
/// and a representative [`LayerHallNumber`] (1..=116) along with the
/// `UnimodularTransformation` that maps the input primitive layer cell to
/// the database canonical for that Hall.
///
/// **v1 scope.** [`LayerSpaceGroup::new`] performs an identity-basis match
/// only -- it does not yet search over `prim_trans_mat` candidates the way
/// `SpaceGroup::new` does. This is sufficient for the M3 round-trip tests
/// (which feed inputs already in canonical layer basis) and for callers who
/// pre-orient their `LayerCell`. Generalised matching with axis-swap
/// candidates and an arithmetic-class pre-filter is deferred to M5.
#[derive(Debug, Clone, Serialize)]
pub struct LayerSpaceGroup {
    pub number: LayerNumber,
    pub hall_number: LayerHallNumber,
    /// Transformation from the input primitive basis to the canonical
    /// primitive basis of `hall_number`.
    pub transformation: UnimodularTransformation,
}

impl LayerSpaceGroup {
    /// Identify a layer group from primitive layer-cell operations.
    ///
    /// Iterates the Hall numbers in `setting` and returns the first one
    /// whose database generators match `prim_layer_operations` modulo an
    /// origin shift. See the type-level docs for the v1 scope (no basis
    /// search).
    pub fn new(
        prim_layer_operations: &Operations,
        setting: LayerSetting,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        let input_order = prim_layer_operations.len();
        for hall_number in setting.hall_numbers() {
            let lh_symbol = LayerHallSymbol::from_hall_number(hall_number)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;

            // Skip groups of the wrong order: without the bulk pipeline's
            // arithmetic-class pre-filter, smaller groups (e.g. LG 1 = p1)
            // would otherwise spuriously match larger inputs since their
            // generators are a strict subset.
            let db_prim_operations = lh_symbol.primitive_traverse();
            if db_prim_operations.len() != input_order {
                continue;
            }

            let db_prim_generators = lh_symbol.primitive_generators();
            if let Some(origin_shift) = match_origin_shift(
                prim_layer_operations,
                &UnimodularLinear::identity(),
                &db_prim_generators,
                epsilon,
            ) {
                let entry = layer_hall_symbol_entry(hall_number)
                    .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
                debug!(
                    "Matched layer Hall number {} (LG {})",
                    hall_number, entry.number
                );
                return Ok(Self {
                    number: entry.number,
                    hall_number,
                    transformation: UnimodularTransformation::new(
                        UnimodularLinear::identity(),
                        origin_shift,
                    ),
                });
            }
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
    /// `LayerSpaceGroup::new`; the returned `(number, hall_number)` must
    /// match the database entry, and the identified transformation must
    /// be unimodular and recover the same operations modulo translation.
    #[test]
    fn test_round_trip_all_layer_hall_numbers() {
        for entry in iter_layer_hall_symbol_entry() {
            let lh_symbol = LayerHallSymbol::from_hall_number(entry.hall_number).unwrap();
            let prim_operations = lh_symbol.primitive_traverse();

            let identified = LayerSpaceGroup::new(
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
                LayerSpaceGroup::new(&prim_operations, LayerSetting::Standard, 1e-8).unwrap();
            assert_eq!(identified.hall_number, hall_number);
        }
    }
}
