use serde::Serialize;

use super::layer_arithmetic_crystal_class::LayerArithmeticNumber;
use super::layer_centering::LayerCentering;

/// Layer-group number (1 - 80), per Fu et al. 2024 Table 5.
pub type LayerNumber = i32;
/// Sequential number for layer-group Hall settings (1 - 116). moyo mints its
/// own numbering since the paper does not number Table 5 rows explicitly; the
/// ordering matches spglib's `database/layer_spg.csv` row order.
pub type LayerHallNumber = i32;

/// One row of the layer-group Hall symbol database.
///
/// Schema mirrors [`super::HallSymbolEntry`] but with [`LayerCentering`]
/// (only `P`/`C` exist for layer groups, paper §3.3) and the layer
/// arithmetic crystal class (1-43, paper Table 4) instead of the bulk one.
#[derive(Debug, Clone, Serialize)]
pub struct LayerHallSymbolEntry {
    /// Sequential layer-group Hall number (1 - 116).
    pub hall_number: LayerHallNumber,
    /// Layer-group number (1 - 80).
    pub number: LayerNumber,
    /// Layer arithmetic crystal class number (1 - 43).
    pub arithmetic_number: LayerArithmeticNumber,
    /// Setting code (paper Table 5 axis/origin labels: "", "a", "b", "b-ac",
    /// "c", "c1", "c2", "c3", "1", "2"). Empty when the LG has a single
    /// setting.
    pub setting: &'static str,
    /// Hall symbol with lowercase `p`/`c` lattice prefix (layer convention).
    pub hall_symbol: &'static str,
    /// Hermann-Mauguin symbol, short notation.
    pub hm_short: &'static str,
    /// Hermann-Mauguin symbol, full notation.
    pub hm_full: &'static str,
    /// Layer centering.
    pub centering: LayerCentering,
}

impl LayerHallSymbolEntry {
    #[allow(clippy::too_many_arguments)]
    const fn new(
        hall_number: LayerHallNumber,
        number: LayerNumber,
        arithmetic_number: LayerArithmeticNumber,
        setting: &'static str,
        hall_symbol: &'static str,
        hm_short: &'static str,
        hm_full: &'static str,
        centering: LayerCentering,
    ) -> Self {
        Self {
            hall_number,
            number,
            arithmetic_number,
            setting,
            hall_symbol,
            hm_short,
            hm_full,
            centering,
        }
    }
}

/// Look up a layer-group Hall entry by `hall_number` (1-based).
pub fn layer_hall_symbol_entry(
    hall_number: LayerHallNumber,
) -> Option<&'static LayerHallSymbolEntry> {
    if hall_number < 1 {
        return None;
    }
    LAYER_HALL_SYMBOL_DATABASE.get((hall_number - 1) as usize)
}

/// Iterate all layer-group Hall entries in `hall_number` order.
pub fn iter_layer_hall_symbol_entry() -> impl Iterator<Item = &'static LayerHallSymbolEntry> {
    LAYER_HALL_SYMBOL_DATABASE.iter()
}

// Source: https://github.com/spglib/spglib/blob/develop/database/layer_spg.csv
// (pinned commit 12355c77fb7c505a55f52cae36341d73b781a065). Regenerate with
// `scripts/layer_hall_symbol.py`. The LG -> layer arithmetic crystal class
// (1-43) mapping is hard-coded in the script and matches the smallest-LG
// boundaries declared in `super::layer_arithmetic_crystal_class`.
#[rustfmt::skip]
const LAYER_HALL_SYMBOL_DATABASE: [LayerHallSymbolEntry; 116] = [
    LayerHallSymbolEntry::new(1, 1, 1, "", "p 1", "p 1", "p 1", LayerCentering::P),
    LayerHallSymbolEntry::new(2, 2, 2, "", "-p 1", "p -1", "p -1", LayerCentering::P),
    LayerHallSymbolEntry::new(3, 3, 3, "c", "p 2", "p 1 1 2", "p 1 1 2", LayerCentering::P),
    LayerHallSymbolEntry::new(4, 4, 4, "c", "p -2", "p 1 1 m", "p 1 1 m", LayerCentering::P),
    LayerHallSymbolEntry::new(5, 5, 4, "c1", "p -2a", "p 1 1 a", "p 1 1 a", LayerCentering::P),
    LayerHallSymbolEntry::new(6, 5, 4, "c2", "p -2ab", "p 1 1 n", "p 1 1 n", LayerCentering::P),
    LayerHallSymbolEntry::new(7, 5, 4, "c3", "p -2b", "p 1 1 b", "p 1 1 b", LayerCentering::P),
    LayerHallSymbolEntry::new(8, 6, 5, "c", "-p 2", "p 1 1 2/m", "p 1 1 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(9, 7, 5, "c1", "-p 2a", "p 1 1 2/a", "p 1 1 2/a", LayerCentering::P),
    LayerHallSymbolEntry::new(10, 7, 5, "c2", "-p 2ab", "p 1 1 2/n", "p 1 1 2/n", LayerCentering::P),
    LayerHallSymbolEntry::new(11, 7, 5, "c3", "-p 2b", "p 1 1 2/b", "p 1 1 2/b", LayerCentering::P),
    LayerHallSymbolEntry::new(12, 8, 6, "a", "p 2x", "p 2 1 1", "p 2 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(13, 8, 6, "b", "p 2y", "p 1 2 1", "p 1 2 1", LayerCentering::P),
    LayerHallSymbolEntry::new(14, 9, 6, "a", "p 2xa", "p 2_1 1 1", "p 2_1 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(15, 9, 6, "b", "p 2yb", "p 1 2_1 1", "p 1 2_1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(16, 10, 7, "a", "c 2x", "c 2 1 1", "c 2 1 1", LayerCentering::C),
    LayerHallSymbolEntry::new(17, 10, 7, "b", "c 2y", "c 1 2 1", "c 1 2 1", LayerCentering::C),
    LayerHallSymbolEntry::new(18, 11, 8, "a", "p -2x", "p m 1 1", "p m 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(19, 11, 8, "b", "p -2y", "p 1 m 1", "p 1 m 1", LayerCentering::P),
    LayerHallSymbolEntry::new(20, 12, 8, "a", "p -2xb", "p b 1 1", "p b 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(21, 12, 8, "b", "p -2ya", "p 1 a 1", "p 1 a 1", LayerCentering::P),
    LayerHallSymbolEntry::new(22, 13, 9, "a", "c -2x", "c m 1 1", "c m 1 1", LayerCentering::C),
    LayerHallSymbolEntry::new(23, 13, 9, "b", "c -2y", "c 1 m 1", "c 1 m 1", LayerCentering::C),
    LayerHallSymbolEntry::new(24, 14, 10, "a", "-p 2x", "p 2/m 1 1", "p 2/m 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(25, 14, 10, "b", "-p 2y", "p 1 2/m 1", "p 1 2/m 1", LayerCentering::P),
    LayerHallSymbolEntry::new(26, 15, 10, "a", "-p 2xa", "p 2_1/m 1 1", "p 2_1/m 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(27, 15, 10, "b", "-p 2yb", "p 1 2_1/m 1", "p 1 2_1/m 1", LayerCentering::P),
    LayerHallSymbolEntry::new(28, 16, 10, "a", "-p 2xb", "p 2/b 1 1", "p 2/b 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(29, 16, 10, "b", "-p 2ya", "p 1 2/a 1", "p 1 2/a 1", LayerCentering::P),
    LayerHallSymbolEntry::new(30, 17, 10, "a", "-p 2xab", "p 2_1/b 1 1", "p 2_1/b 1 1", LayerCentering::P),
    LayerHallSymbolEntry::new(31, 17, 10, "b", "-p 2yab", "p 1 2_1/a 1", "p 1 2_1/a 1", LayerCentering::P),
    LayerHallSymbolEntry::new(32, 18, 11, "a", "-c 2x", "c 2/m 1 1", "c 2/m 1 1", LayerCentering::C),
    LayerHallSymbolEntry::new(33, 18, 11, "b", "-c 2y", "c 1 2/m 1", "c 1 2/m 1", LayerCentering::C),
    LayerHallSymbolEntry::new(34, 19, 12, "", "p 2 2", "p 2 2 2", "p 2 2 2", LayerCentering::P),
    LayerHallSymbolEntry::new(35, 20, 12, "", "p 2 2a", "p 2_1 2 2", "p 2_1 2 2", LayerCentering::P),
    LayerHallSymbolEntry::new(36, 20, 12, "b-ac", "p 2 2b", "p 2 2_1 2", "p 2 2_1 2", LayerCentering::P),
    LayerHallSymbolEntry::new(37, 21, 12, "", "p 2 2ab", "p 2_1 2_1 2", "p 2_1 2_1 2", LayerCentering::P),
    LayerHallSymbolEntry::new(38, 22, 13, "", "c 2 2", "c 2 2 2", "c 2 2 2", LayerCentering::C),
    LayerHallSymbolEntry::new(39, 23, 14, "", "p 2 -2", "p m m 2", "p m m 2", LayerCentering::P),
    LayerHallSymbolEntry::new(40, 24, 14, "", "p 2 -2a", "p m a 2", "p m a 2", LayerCentering::P),
    LayerHallSymbolEntry::new(41, 24, 14, "b-ac", "p 2 -2b", "p b m 2", "p b m 2", LayerCentering::P),
    LayerHallSymbolEntry::new(42, 25, 14, "", "p 2 -2ab", "p b a 2", "p b a 2", LayerCentering::P),
    LayerHallSymbolEntry::new(43, 26, 15, "", "c 2 -2", "c m m 2", "c m m 2", LayerCentering::C),
    LayerHallSymbolEntry::new(44, 27, 16, "", "p -2 -2", "p m 2 m", "p m 2 m", LayerCentering::P),
    LayerHallSymbolEntry::new(45, 27, 16, "b-ac", "p -2 2", "p 2 m m", "p 2 m m", LayerCentering::P),
    LayerHallSymbolEntry::new(46, 28, 16, "", "p -2b -2", "p m 2_1 b", "p m 2_1 b", LayerCentering::P),
    LayerHallSymbolEntry::new(47, 28, 16, "b-ac", "p -2a 2a", "p 2_1 m a", "p 2_1 m a", LayerCentering::P),
    LayerHallSymbolEntry::new(48, 29, 16, "", "p -2 -2b", "p b 2_1 m", "p b 2_1 m", LayerCentering::P),
    LayerHallSymbolEntry::new(49, 29, 16, "b-ac", "p -2 2a", "p 2_1 a m", "p 2_1 a m", LayerCentering::P),
    LayerHallSymbolEntry::new(50, 30, 16, "", "p -2b -2b", "p b 2 b", "p b 2 b", LayerCentering::P),
    LayerHallSymbolEntry::new(51, 30, 16, "b-ac", "p -2a 2", "p 2 a a", "p 2 a a", LayerCentering::P),
    LayerHallSymbolEntry::new(52, 31, 16, "", "p -2a -2a", "p m 2 a", "p m 2 a", LayerCentering::P),
    LayerHallSymbolEntry::new(53, 31, 16, "b-ac", "p -2b 2", "p 2 m b", "p 2 m b", LayerCentering::P),
    LayerHallSymbolEntry::new(54, 32, 16, "", "p -2ab -2", "p m 2_1 n", "p m 2_1 n", LayerCentering::P),
    LayerHallSymbolEntry::new(55, 32, 16, "b-ac", "p -2ab 2ab", "p 2_1 m n", "p 2_1 m n", LayerCentering::P),
    LayerHallSymbolEntry::new(56, 33, 16, "", "p -2a -2ab", "p b 2_1 a", "p b 2_1 a", LayerCentering::P),
    LayerHallSymbolEntry::new(57, 33, 16, "b-ac", "p -2b 2a", "p 2_1 a b", "p 2_1 a b", LayerCentering::P),
    LayerHallSymbolEntry::new(58, 34, 16, "", "p -2ab -2ab", "p b 2 n", "p b 2 n", LayerCentering::P),
    LayerHallSymbolEntry::new(59, 34, 16, "b-ac", "p -2ab 2", "p 2 a n", "p 2 a n", LayerCentering::P),
    LayerHallSymbolEntry::new(60, 35, 17, "", "c -2 -2", "c m 2 m", "c m 2 m", LayerCentering::C),
    LayerHallSymbolEntry::new(61, 35, 17, "b-ac", "c -2 2", "c 2 m m", "c 2 m m", LayerCentering::C),
    LayerHallSymbolEntry::new(62, 36, 17, "", "c -2a -2a", "c m 2 e", "c m 2 e", LayerCentering::C),
    LayerHallSymbolEntry::new(63, 36, 17, "b-ac", "c -2a 2", "c 2 m e", "c 2 m e", LayerCentering::C),
    LayerHallSymbolEntry::new(64, 37, 18, "", "-p 2 2", "p m m m", "p 2/m 2/m 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(65, 38, 18, "", "-p 2a 2", "p m a a", "p 2/m 2/a 2/a", LayerCentering::P),
    LayerHallSymbolEntry::new(66, 38, 18, "b-ac", "-p 2b 2b", "p b m b", "p 2/b 2/m 2/b", LayerCentering::P),
    LayerHallSymbolEntry::new(67, 39, 18, "", "-p 2ab 2b", "p b a n", "p 2/b 2/a 2/n", LayerCentering::P),
    LayerHallSymbolEntry::new(68, 40, 18, "", "-p 2 2a", "p m a m", "p 2_1/m 2/a 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(69, 40, 18, "b-ac", "-p 2 2b", "p b m m", "p 2/b 2_1/m 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(70, 41, 18, "", "-p 2a 2a", "p m m a", "p 2_1/m 2/m 2/a", LayerCentering::P),
    LayerHallSymbolEntry::new(71, 41, 18, "b-ac", "-p 2b 2", "p m m b", "p 2/m 2_1/m 2/b", LayerCentering::P),
    LayerHallSymbolEntry::new(72, 42, 18, "", "-p 2ab 2", "p m a n", "p 2/m 2_1/a 2/n", LayerCentering::P),
    LayerHallSymbolEntry::new(73, 42, 18, "b-ac", "-p 2ab 2ab", "p b m n", "p 2_1/b 2/m 2/n", LayerCentering::P),
    LayerHallSymbolEntry::new(74, 43, 18, "", "-p 2a 2b", "p b a a", "p 2/b 2_1/a 2/a", LayerCentering::P),
    LayerHallSymbolEntry::new(75, 43, 18, "b-ac", "-p 2b 2ab", "p b a b", "p 2_1/b 2/a 2/b", LayerCentering::P),
    LayerHallSymbolEntry::new(76, 44, 18, "", "-p 2 2ab", "p b a m", "p 2_1/b 2_1/a 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(77, 45, 18, "", "-p 2a 2ab", "p b m a", "p 2_1/b 2_1/m 2/a", LayerCentering::P),
    LayerHallSymbolEntry::new(78, 45, 18, "b-ac", "-p 2b 2a", "p m a b", "p 2_1/m 2_1/a 2/b", LayerCentering::P),
    LayerHallSymbolEntry::new(79, 46, 18, "", "-p 2ab 2a", "p m m n", "p 2_1/m 2_1/m 2/n", LayerCentering::P),
    LayerHallSymbolEntry::new(80, 47, 19, "", "-c 2 2", "c m m m", "c 2/m 2/m 2/m", LayerCentering::C),
    LayerHallSymbolEntry::new(81, 48, 19, "", "-c 2a 2", "c m m e", "c 2/m 2/m 2/e", LayerCentering::C),
    LayerHallSymbolEntry::new(82, 49, 20, "", "p 4", "p 4", "p 4", LayerCentering::P),
    LayerHallSymbolEntry::new(83, 50, 21, "", "p -4", "p -4", "p -4", LayerCentering::P),
    LayerHallSymbolEntry::new(84, 51, 22, "", "-p 4", "p 4/m", "p 4/m", LayerCentering::P),
    LayerHallSymbolEntry::new(85, 52, 22, "1", "p 4 -1ab", "p 4/n", "p 4/n", LayerCentering::P),
    LayerHallSymbolEntry::new(86, 52, 22, "2", "-p 4a", "p 4/n", "p 4/n", LayerCentering::P),
    LayerHallSymbolEntry::new(87, 53, 23, "", "p 4 2", "p 4 2 2", "p 4 2 2", LayerCentering::P),
    LayerHallSymbolEntry::new(88, 54, 23, "", "p 4 2ab", "p 4 2_1 2", "p 4 2_1 2", LayerCentering::P),
    LayerHallSymbolEntry::new(89, 55, 24, "", "p 4 -2", "p 4 m m", "p 4 m m", LayerCentering::P),
    LayerHallSymbolEntry::new(90, 56, 24, "", "p 4 -2ab", "p 4 b m", "p 4 b m", LayerCentering::P),
    LayerHallSymbolEntry::new(91, 57, 25, "", "p -4 2", "p -4 2 m", "p -4 2 m", LayerCentering::P),
    LayerHallSymbolEntry::new(92, 58, 25, "", "p -4 2ab", "p -4 2_1 m", "p -4 2_1 m", LayerCentering::P),
    LayerHallSymbolEntry::new(93, 59, 26, "", "p -4 -2", "p -4 m 2", "p -4 m 2", LayerCentering::P),
    LayerHallSymbolEntry::new(94, 60, 26, "", "p -4 -2ab", "p -4 b 2", "p -4 b 2", LayerCentering::P),
    LayerHallSymbolEntry::new(95, 61, 27, "", "-p 4 2", "p 4/m m m", "p 4/m 2/m 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(96, 62, 27, "1", "p 4 2 -1ab", "p 4/n b m", "p 4/n 2/b 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(97, 62, 27, "2", "-p 4a 2b", "p 4/n b m", "p 4/n 2/b 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(98, 63, 27, "", "-p 4 2ab", "p 4/m b m", "p 4/m 2_1/b m", LayerCentering::P),
    LayerHallSymbolEntry::new(99, 64, 27, "1", "p 4 2ab -1ab", "p 4/n m m", "p 4/n 2_1/m m", LayerCentering::P),
    LayerHallSymbolEntry::new(100, 64, 27, "2", "-p 4a 2a", "p 4/n m m", "p 4/n 2_1/m m", LayerCentering::P),
    LayerHallSymbolEntry::new(101, 65, 28, "", "p 3", "p 3", "p 3", LayerCentering::P),
    LayerHallSymbolEntry::new(102, 66, 29, "", "-p 3", "p -3", "p -3", LayerCentering::P),
    LayerHallSymbolEntry::new(103, 67, 30, "", "p 3 2", "p 3 1 2", "p 3 1 2", LayerCentering::P),
    LayerHallSymbolEntry::new(104, 68, 31, "", "p 3 2=", "p 3 2 1", "p 3 2 1", LayerCentering::P),
    LayerHallSymbolEntry::new(105, 69, 32, "", "p 3 -2=", "p 3 m 1", "p 3 m 1", LayerCentering::P),
    LayerHallSymbolEntry::new(106, 70, 33, "", "p 3 -2", "p 3 1 m", "p 3 1 m", LayerCentering::P),
    LayerHallSymbolEntry::new(107, 71, 34, "", "-p 3 2", "p -3 1 m", "p -3 1 2/m", LayerCentering::P),
    LayerHallSymbolEntry::new(108, 72, 35, "", "-p 3 2=", "p -3 m 1", "p -3 2/m 1", LayerCentering::P),
    LayerHallSymbolEntry::new(109, 73, 36, "", "p 6", "p 6", "p 6", LayerCentering::P),
    LayerHallSymbolEntry::new(110, 74, 37, "", "p -6", "p -6", "p -6", LayerCentering::P),
    LayerHallSymbolEntry::new(111, 75, 38, "", "-p 6", "p 6/m", "p 6/m", LayerCentering::P),
    LayerHallSymbolEntry::new(112, 76, 39, "", "p 6 2", "p 6 2 2", "p 6 2 2", LayerCentering::P),
    LayerHallSymbolEntry::new(113, 77, 40, "", "p 6 -2", "p 6 m m", "p 6 m m", LayerCentering::P),
    LayerHallSymbolEntry::new(114, 78, 41, "", "p -6 2", "p -6 m 2", "p -6 m 2", LayerCentering::P),
    LayerHallSymbolEntry::new(115, 79, 42, "", "p -6 -2", "p -6 2 m", "p -6 2 m", LayerCentering::P),
    LayerHallSymbolEntry::new(116, 80, 43, "", "-p 6 2", "p 6/m m m", "p 6/m m m", LayerCentering::P),
];

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn test_database_size_and_numbering() {
        assert_eq!(LAYER_HALL_SYMBOL_DATABASE.len(), 116);
        for (i, entry) in iter_layer_hall_symbol_entry().enumerate() {
            assert_eq!(entry.hall_number as usize, i + 1);
        }
    }

    #[test]
    fn test_lookup_bounds() {
        assert!(layer_hall_symbol_entry(0).is_none());
        assert!(layer_hall_symbol_entry(117).is_none());
        assert!(layer_hall_symbol_entry(1).is_some());
        assert!(layer_hall_symbol_entry(116).is_some());
    }

    #[test]
    fn test_every_lg_covered() {
        let mut covered: HashSet<LayerNumber> = HashSet::new();
        for entry in iter_layer_hall_symbol_entry() {
            covered.insert(entry.number);
        }
        for lg in 1..=80 {
            assert!(covered.contains(&lg), "LG {} missing", lg);
        }
    }

    #[test]
    fn test_centering_matches_hall_prefix() {
        for entry in iter_layer_hall_symbol_entry() {
            let head = entry.hall_symbol.trim_start_matches('-').trim_start();
            let expected = if head.starts_with("c ") {
                LayerCentering::C
            } else if head.starts_with("p ") {
                LayerCentering::P
            } else {
                panic!("unexpected layer Hall prefix: {}", entry.hall_symbol);
            };
            assert_eq!(
                entry.centering, expected,
                "centering mismatch for hall_number {}",
                entry.hall_number
            );
        }
    }
}
