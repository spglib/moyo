use std::collections::HashMap;

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

// In-plane normalizer corrections (besides identity) that preserve the
// layer block form `W_i3 = W_3i = 0, W_33 = ±1` (paper Fu et al. 2024
// eq. 4).
//
// * `b a -c` is the in-plane `a` ↔ `b` swap paired with a `c` sign-flip;
//   it covers monoclinic-rectangular `:b` ↔ `:a` and orthorhombic
//   `b-ac` ↔ `abc` alignment when the input lies in a non-canonical
//   Hall basis.
// * `-a b -c` covers the trigonal / hexagonal axis-orientation
//   ambiguity: a hexagonal in-plane lattice has two inequivalent
//   choices of (a, b) that share the same lengths and 120 deg angle,
//   and conjugation by `diag(-1, +1, -1)` swaps the two 3-fold
//   conventions (`[[0, -1, 0], [1, -1, 0], [0, 0, 1]]` vs
//   `[[0, 1, 0], [-1, -1, 0], [0, 0, 1]]`). Without this correction LG
//   71 / 72 (and their hexagonal supergroups) silently fail to
//   identify on inputs whose primitive basis matches the "other" hex
//   convention.
//
// For `c`-centered layer Bravais classes (`oc`) the primitive cell of
// a C-centered conventional admits multiple equivalent choices that
// differ by an in-plane shear of infinite order in GL(2, Z); a
// finite correction-matrix list cannot cover this. The proper fix
// mirrors the bulk `correction_transformation_matrices` flow --
// take the bulk monoclinic / orthorhombic conventional `convs` and
// conjugate them by the layer centering transform -- and is left as
// future work. JARVIS bench failures dominated by SG 12 / 21 / 65
// (LG 18, LG 22, LG 26, ...) hit this gap.
const LAYER_CORRECTION_MATRICES: [UnimodularLinear; 3] = [
    UnimodularLinear::new(
        1, 0, 0, //
        0, 1, 0, //
        0, 0, 1, //
    ),
    UnimodularLinear::new(
        0, 1, 0, //
        1, 0, 0, //
        0, 0, -1, //
    ),
    UnimodularLinear::new(
        -1, 0, 0, //
        0, 1, 0, //
        0, 0, -1, //
    ),
];

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
        // Mirror `SpaceGroup::new` but consume only the geometric class:
        // the bulk arithmetic class normalises to a 3D-canonical orientation
        // that re-labels the layer's c-axis, so its `prim_trans_mat` is not
        // re-used here.
        let prim_rotations = project_rotations(prim_layer_operations);
        let point_group = PointGroup::new(&prim_rotations)?;
        let geometric_crystal_class = arithmetic_crystal_class_entry(point_group.arithmetic_number)
            .unwrap()
            .geometric_crystal_class;
        debug!(
            "Layer point group: geometric crystal class {:?}",
            geometric_crystal_class
        );

        for hall_number in setting.hall_numbers() {
            let entry = layer_hall_symbol_entry(hall_number).unwrap();
            let layer_arith =
                layer_arithmetic_crystal_class_entry(entry.arithmetic_number).unwrap();
            if layer_arith.geometric_crystal_class != geometric_crystal_class {
                continue;
            }

            let lh_symbol = LayerHallSymbol::new(entry.hall_symbol)
                .ok_or(MoyoError::LayerGroupTypeIdentificationError)?;
            let db_prim_generators = lh_symbol.primitive_generators();

            for trans_mat in &LAYER_CORRECTION_MATRICES {
                if let Some(mut origin_shift) = match_origin_shift(
                    prim_layer_operations,
                    trans_mat,
                    &db_prim_generators,
                    epsilon,
                ) {
                    // The c-axis is aperiodic, but `match_origin_shift` solves
                    // its modular system in all three components. For LGs
                    // with c-flipping generators (`W[2,2] = -1`) the c-row
                    // reduces to `-2 s_z = b_z (mod 1)` and admits two valid
                    // branches differing by 1/2; the modular solver returns
                    // either, but only one places the layer's special-z
                    // atoms onto the LG Wyckoff database orbits (stored at
                    // `z = 0` -- there is no `z = 1/2` counterpart for
                    // layer groups). The layer search step now produces
                    // operations whose `t_z` is exact (no mod-1 reduction),
                    // so we can recover `s_z` directly from the c-row of any
                    // c-flipping generator without going through `% 1`.
                    //
                    // For any layer-block W with `W[2,2] = -1`, `(W - E)`
                    // has c-row `(0, 0, -2)`, so `-2 s_z = t_db_z - t_target_z`
                    // and `s_z = (t_target_z - t_db_z) / 2` exactly.
                    // `LAYER_CORRECTION_MATRICES` are block-diagonal in c
                    // (`trans_mat[2,0] = trans_mat[2,1] = 0`,
                    // `trans_mat[2,2] in {+1, -1}`) so the c-component of
                    // `trans_mat * s` is `trans_mat[2,2] * s_z`.
                    if let Some(exact_s_z) = layer_exact_s_z(
                        prim_layer_operations,
                        trans_mat,
                        &db_prim_generators,
                        epsilon,
                    ) {
                        origin_shift[2] = trans_mat[(2, 2)] as f64 * exact_s_z;
                    }
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
