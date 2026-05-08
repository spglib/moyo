use super::super::classification::GeometricCrystalClass;
use super::layer_classification::{LayerBravaisClass, LayerLatticeSystem};

pub type LayerArithmeticNumber = i32;

#[derive(Debug, Clone)]
pub struct LayerArithmeticCrystalClassEntry {
    /// Sequential number for layer arithmetic crystal classes (1-based, ordered as in paper Table 4).
    pub arithmetic_number: LayerArithmeticNumber,
    /// Symbol for the arithmetic crystal class (e.g. "p1", "c2/m11").
    pub symbol: &'static str,
    /// Geometric crystal class (subset of the 32 used for space groups; cubic classes do not occur).
    pub geometric_crystal_class: GeometricCrystalClass,
    /// Bravais class for the layer group's 2D lattice.
    pub layer_bravais_class: LayerBravaisClass,
}

impl LayerArithmeticCrystalClassEntry {
    const fn new(
        arithmetic_number: LayerArithmeticNumber,
        symbol: &'static str,
        geometric_crystal_class: GeometricCrystalClass,
        layer_bravais_class: LayerBravaisClass,
    ) -> Self {
        Self {
            arithmetic_number,
            symbol,
            geometric_crystal_class,
            layer_bravais_class,
        }
    }

    pub fn layer_lattice_system(&self) -> LayerLatticeSystem {
        LayerLatticeSystem::from_layer_bravais_class(self.layer_bravais_class)
    }
}

pub fn layer_arithmetic_crystal_class_entry(
    arithmetic_number: LayerArithmeticNumber,
) -> Option<LayerArithmeticCrystalClassEntry> {
    if arithmetic_number < 1 {
        return None;
    }
    LAYER_ARITHMETIC_CRYSTAL_CLASS_DATABASE
        .get((arithmetic_number - 1) as usize)
        .cloned()
}

pub fn iter_layer_arithmetic_crystal_entry()
-> impl Iterator<Item = &'static LayerArithmeticCrystalClassEntry> {
    LAYER_ARITHMETIC_CRYSTAL_CLASS_DATABASE.iter()
}

// Layer arithmetic-geometric crystal classes from Fu et al. 2024, Table 4.
// Each row of Table 4 is one arithmetic class (orientation-distinct rows like
// "(m2m)" vs "mm2" are listed separately, mirroring how the existing 3D
// `ArithmeticCrystalClassEntry` distinguishes e.g. `mm2C` and `2mmC`).
// The integer in parentheses after each symbol below is the smallest layer-group
// number (LG number) belonging to the class, as listed in Table 4.
const LAYER_ARITHMETIC_CRYSTAL_CLASS_DATABASE: [LayerArithmeticCrystalClassEntry; 43] = [
    // Triclinic / oblique
    LayerArithmeticCrystalClassEntry::new(
        1,
        "p1",
        GeometricCrystalClass::C1,
        LayerBravaisClass::mp,
    ), // LG 1
    LayerArithmeticCrystalClassEntry::new(
        2,
        "p-1",
        GeometricCrystalClass::Ci,
        LayerBravaisClass::mp,
    ), // LG 2
    // Monoclinic / oblique
    LayerArithmeticCrystalClassEntry::new(
        3,
        "p112",
        GeometricCrystalClass::C2,
        LayerBravaisClass::mp,
    ), // LG 3
    LayerArithmeticCrystalClassEntry::new(
        4,
        "p11m",
        GeometricCrystalClass::C1h,
        LayerBravaisClass::mp,
    ), // LG 4
    LayerArithmeticCrystalClassEntry::new(
        5,
        "p112/m",
        GeometricCrystalClass::C2h,
        LayerBravaisClass::mp,
    ), // LG 6
    // Monoclinic / rectangular
    LayerArithmeticCrystalClassEntry::new(
        6,
        "p211",
        GeometricCrystalClass::C2,
        LayerBravaisClass::op,
    ), // LG 8
    LayerArithmeticCrystalClassEntry::new(
        7,
        "c211",
        GeometricCrystalClass::C2,
        LayerBravaisClass::oc,
    ), // LG 10
    LayerArithmeticCrystalClassEntry::new(
        8,
        "pm11",
        GeometricCrystalClass::C1h,
        LayerBravaisClass::op,
    ), // LG 11
    LayerArithmeticCrystalClassEntry::new(
        9,
        "cm11",
        GeometricCrystalClass::C1h,
        LayerBravaisClass::oc,
    ), // LG 13
    LayerArithmeticCrystalClassEntry::new(
        10,
        "p2/m11",
        GeometricCrystalClass::C2h,
        LayerBravaisClass::op,
    ), // LG 14
    LayerArithmeticCrystalClassEntry::new(
        11,
        "c2/m11",
        GeometricCrystalClass::C2h,
        LayerBravaisClass::oc,
    ), // LG 18
    // Orthorhombic / rectangular
    LayerArithmeticCrystalClassEntry::new(
        12,
        "p222",
        GeometricCrystalClass::D2,
        LayerBravaisClass::op,
    ), // LG 19
    LayerArithmeticCrystalClassEntry::new(
        13,
        "c222",
        GeometricCrystalClass::D2,
        LayerBravaisClass::oc,
    ), // LG 22
    LayerArithmeticCrystalClassEntry::new(
        14,
        "pmm2",
        GeometricCrystalClass::C2v,
        LayerBravaisClass::op,
    ), // LG 23
    LayerArithmeticCrystalClassEntry::new(
        15,
        "cmm2",
        GeometricCrystalClass::C2v,
        LayerBravaisClass::oc,
    ), // LG 26
    LayerArithmeticCrystalClassEntry::new(
        16,
        "pm2m",
        GeometricCrystalClass::C2v,
        LayerBravaisClass::op,
    ), // LG 27
    LayerArithmeticCrystalClassEntry::new(
        17,
        "cm2m",
        GeometricCrystalClass::C2v,
        LayerBravaisClass::oc,
    ), // LG 35
    LayerArithmeticCrystalClassEntry::new(
        18,
        "pmmm",
        GeometricCrystalClass::D2h,
        LayerBravaisClass::op,
    ), // LG 37
    LayerArithmeticCrystalClassEntry::new(
        19,
        "cmmm",
        GeometricCrystalClass::D2h,
        LayerBravaisClass::oc,
    ), // LG 47
    // Tetragonal / square
    LayerArithmeticCrystalClassEntry::new(
        20,
        "p4",
        GeometricCrystalClass::C4,
        LayerBravaisClass::tp,
    ), // LG 49
    LayerArithmeticCrystalClassEntry::new(
        21,
        "p-4",
        GeometricCrystalClass::S4,
        LayerBravaisClass::tp,
    ), // LG 50
    LayerArithmeticCrystalClassEntry::new(
        22,
        "p4/m",
        GeometricCrystalClass::C4h,
        LayerBravaisClass::tp,
    ), // LG 51
    LayerArithmeticCrystalClassEntry::new(
        23,
        "p422",
        GeometricCrystalClass::D4,
        LayerBravaisClass::tp,
    ), // LG 53
    LayerArithmeticCrystalClassEntry::new(
        24,
        "p4mm",
        GeometricCrystalClass::C4v,
        LayerBravaisClass::tp,
    ), // LG 55
    LayerArithmeticCrystalClassEntry::new(
        25,
        "p-42m",
        GeometricCrystalClass::D2d,
        LayerBravaisClass::tp,
    ), // LG 57
    LayerArithmeticCrystalClassEntry::new(
        26,
        "p-4m2",
        GeometricCrystalClass::D2d,
        LayerBravaisClass::tp,
    ), // LG 59
    LayerArithmeticCrystalClassEntry::new(
        27,
        "p4/mmm",
        GeometricCrystalClass::D4h,
        LayerBravaisClass::tp,
    ), // LG 61
    // Trigonal / hexagonal
    LayerArithmeticCrystalClassEntry::new(
        28,
        "p3",
        GeometricCrystalClass::C3,
        LayerBravaisClass::hp,
    ), // LG 65
    LayerArithmeticCrystalClassEntry::new(
        29,
        "p-3",
        GeometricCrystalClass::C3i,
        LayerBravaisClass::hp,
    ), // LG 66
    LayerArithmeticCrystalClassEntry::new(
        30,
        "p312",
        GeometricCrystalClass::D3,
        LayerBravaisClass::hp,
    ), // LG 67
    LayerArithmeticCrystalClassEntry::new(
        31,
        "p321",
        GeometricCrystalClass::D3,
        LayerBravaisClass::hp,
    ), // LG 68
    LayerArithmeticCrystalClassEntry::new(
        32,
        "p3m1",
        GeometricCrystalClass::C3v,
        LayerBravaisClass::hp,
    ), // LG 69
    LayerArithmeticCrystalClassEntry::new(
        33,
        "p31m",
        GeometricCrystalClass::C3v,
        LayerBravaisClass::hp,
    ), // LG 70
    LayerArithmeticCrystalClassEntry::new(
        34,
        "p-31m",
        GeometricCrystalClass::D3d,
        LayerBravaisClass::hp,
    ), // LG 71
    LayerArithmeticCrystalClassEntry::new(
        35,
        "p-3m1",
        GeometricCrystalClass::D3d,
        LayerBravaisClass::hp,
    ), // LG 72
    // Hexagonal / hexagonal
    LayerArithmeticCrystalClassEntry::new(
        36,
        "p6",
        GeometricCrystalClass::C6,
        LayerBravaisClass::hp,
    ), // LG 73
    LayerArithmeticCrystalClassEntry::new(
        37,
        "p-6",
        GeometricCrystalClass::C3h,
        LayerBravaisClass::hp,
    ), // LG 74
    LayerArithmeticCrystalClassEntry::new(
        38,
        "p6/m",
        GeometricCrystalClass::C6h,
        LayerBravaisClass::hp,
    ), // LG 75
    LayerArithmeticCrystalClassEntry::new(
        39,
        "p622",
        GeometricCrystalClass::D6,
        LayerBravaisClass::hp,
    ), // LG 76
    LayerArithmeticCrystalClassEntry::new(
        40,
        "p6mm",
        GeometricCrystalClass::C6v,
        LayerBravaisClass::hp,
    ), // LG 77
    LayerArithmeticCrystalClassEntry::new(
        41,
        "p-6m2",
        GeometricCrystalClass::D3h,
        LayerBravaisClass::hp,
    ), // LG 78
    LayerArithmeticCrystalClassEntry::new(
        42,
        "p-62m",
        GeometricCrystalClass::D3h,
        LayerBravaisClass::hp,
    ), // LG 79
    LayerArithmeticCrystalClassEntry::new(
        43,
        "p6/mmm",
        GeometricCrystalClass::D6h,
        LayerBravaisClass::hp,
    ), // LG 80
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_database_size_and_numbering() {
        assert_eq!(LAYER_ARITHMETIC_CRYSTAL_CLASS_DATABASE.len(), 43);
        for (i, entry) in iter_layer_arithmetic_crystal_entry().enumerate() {
            assert_eq!(entry.arithmetic_number as usize, i + 1);
        }
    }

    #[test]
    fn test_lookup() {
        let entry = layer_arithmetic_crystal_class_entry(1).unwrap();
        assert_eq!(entry.symbol, "p1");
        assert_eq!(entry.geometric_crystal_class, GeometricCrystalClass::C1);
        assert_eq!(entry.layer_bravais_class, LayerBravaisClass::mp);
        assert_eq!(entry.layer_lattice_system(), LayerLatticeSystem::Oblique);

        let entry = layer_arithmetic_crystal_class_entry(43).unwrap();
        assert_eq!(entry.symbol, "p6/mmm");
        assert_eq!(entry.geometric_crystal_class, GeometricCrystalClass::D6h);
        assert_eq!(entry.layer_bravais_class, LayerBravaisClass::hp);

        assert!(layer_arithmetic_crystal_class_entry(0).is_none());
        assert!(layer_arithmetic_crystal_class_entry(44).is_none());
    }

    #[test]
    fn test_no_cubic() {
        for entry in iter_layer_arithmetic_crystal_entry() {
            assert!(!matches!(
                entry.geometric_crystal_class,
                GeometricCrystalClass::T
                    | GeometricCrystalClass::Th
                    | GeometricCrystalClass::O
                    | GeometricCrystalClass::Td
                    | GeometricCrystalClass::Oh
            ));
        }
    }
}
