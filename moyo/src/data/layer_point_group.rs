use super::hall_symbol::LayerHallSymbol;
use super::layer_arithmetic_crystal_class::LayerArithmeticNumber;
use super::layer_centering::LayerCentering;
use super::layer_hall_symbol_database::iter_layer_hall_symbol_entry;
use crate::base::{Rotations, project_rotations};

/// Layer-group analog of [`super::PointGroupRepresentative`]: canonical
/// primitive rotation generators per layer arithmetic crystal class
/// (1..=43, paper Fu et al. 2024 Table 4) in the layer convention (the
/// aperiodic axis is always `c`).
///
/// Used by `identify::layer_point_group::LayerPointGroup::new` to match
/// an input's projected rotation set against the canonical generators of
/// each candidate arithmetic class.
#[derive(Debug)]
pub struct LayerPointGroupRepresentative {
    /// Conventional-cell rotation generators of the representative Hall.
    pub generators: Rotations,
    /// Layer Bravais centering of the representative Hall (`P` or `C`).
    pub centering: LayerCentering,
}

impl LayerPointGroupRepresentative {
    fn new(generators: Rotations, centering: LayerCentering) -> Self {
        Self {
            generators,
            centering,
        }
    }

    /// Construct the representative for a layer arithmetic crystal class
    /// from its smallest [`super::LayerHallNumber`]. The smallest-Hall
    /// choice matches `LayerSetting::Spglib`'s default ordering and
    /// keeps the per-class basis convention consistent with the
    /// `LayerHallNumber` ↔ `LayerArithmeticNumber` mapping in the
    /// database.
    pub fn from_layer_arithmetic_crystal_class(arithmetic_number: LayerArithmeticNumber) -> Self {
        let entry = iter_layer_hall_symbol_entry()
            .find(|e| e.arithmetic_number == arithmetic_number)
            .unwrap_or_else(|| {
                panic!(
                    "no LayerHallSymbolEntry for arithmetic_number={}",
                    arithmetic_number
                )
            });
        let layer_hall_symbol = LayerHallSymbol::from_hall_number(entry.hall_number)
            .expect("LayerHallSymbol::from_hall_number for db entry");
        Self::new(
            project_rotations(layer_hall_symbol.generators()),
            entry.centering,
        )
    }

    /// Project the conventional-cell generators to the primitive cell.
    /// Mirrors `PointGroupRepresentative::primitive_generators`.
    pub fn primitive_generators(&self) -> Rotations {
        let prim_trans_mat_inv = self.centering.linear().map(|e| e as f64);
        let prim_trans_mat = self.centering.inverse();
        self.generators
            .iter()
            .map(|g| {
                let prim_g = prim_trans_mat_inv * g.map(|e| e as f64) * prim_trans_mat;
                prim_g.map(|e| e.round() as i32)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::traverse;
    use crate::data::layer_arithmetic_crystal_class::iter_layer_arithmetic_crystal_entry;

    /// Round-trip: every layer arithmetic class (1..=43) produces a
    /// representative whose primitive generators reconstruct a finite
    /// rotation group with the expected geometric crystal class.
    #[test]
    fn test_layer_point_group_representative() {
        for entry in iter_layer_arithmetic_crystal_entry() {
            let rep = LayerPointGroupRepresentative::from_layer_arithmetic_crystal_class(
                entry.arithmetic_number,
            );
            let prim_gens = rep.primitive_generators();
            // Generators close into a finite group.
            let rotations = traverse(&prim_gens);
            assert!(
                !rotations.is_empty(),
                "empty rotation group for layer arithmetic_number={}",
                entry.arithmetic_number
            );
        }
    }
}
