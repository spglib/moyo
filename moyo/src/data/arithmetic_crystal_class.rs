use super::classification::{BravaisClass, GeometricCrystalClass, LatticeSystem};

pub type ArithmeticNumber = i32;

#[derive(Debug, Clone)]
pub struct ArithmeticCrystalClassEntry {
    /// Number for arithmetic crystal classes (1 - 73)
    pub arithmetic_number: ArithmeticNumber,
    /// Symbol for arithmetic crystal class
    pub symbol: &'static str,
    /// Geometric crystal class
    pub geometric_crystal_class: GeometricCrystalClass,
    /// Bravais class
    pub bravais_class: BravaisClass,
}

impl ArithmeticCrystalClassEntry {
    const fn new(
        arithmetic_number: ArithmeticNumber,
        symbol: &'static str,
        geometric_crystal_class: GeometricCrystalClass,
        bravais_class: BravaisClass,
    ) -> Self {
        Self {
            arithmetic_number,
            symbol,
            geometric_crystal_class,
            bravais_class,
        }
    }

    pub fn lattice_system(&self) -> LatticeSystem {
        LatticeSystem::from_bravais_class(self.bravais_class)
    }
}

pub fn arithmetic_crystal_class_entry(
    arithmetic_number: ArithmeticNumber,
) -> Option<ArithmeticCrystalClassEntry> {
    ARITHMETIC_CRYSTAL_CLASS_DATABASE
        .get(arithmetic_number as usize - 1)
        .cloned()
}

pub fn iter_arithmetic_crystal_entry() -> impl Iterator<Item = &'static ArithmeticCrystalClassEntry>
{
    ARITHMETIC_CRYSTAL_CLASS_DATABASE.iter()
}

// Ordered the same as https://dictionary.iucr.org/Arithmetic_crystal_class
const ARITHMETIC_CRYSTAL_CLASS_DATABASE: [ArithmeticCrystalClassEntry; 73] = [
    // Crystal system: triclinic
    ArithmeticCrystalClassEntry::new(1, "1P", GeometricCrystalClass::C1, BravaisClass::aP),
    ArithmeticCrystalClassEntry::new(2, "-1P", GeometricCrystalClass::Ci, BravaisClass::aP),
    // Crystal system: monoclinic
    ArithmeticCrystalClassEntry::new(3, "2P", GeometricCrystalClass::C2, BravaisClass::mP),
    ArithmeticCrystalClassEntry::new(4, "2C", GeometricCrystalClass::C2, BravaisClass::mC),
    ArithmeticCrystalClassEntry::new(5, "mP", GeometricCrystalClass::C1h, BravaisClass::mP),
    ArithmeticCrystalClassEntry::new(6, "mC", GeometricCrystalClass::C1h, BravaisClass::mC),
    ArithmeticCrystalClassEntry::new(7, "2/mP", GeometricCrystalClass::C2h, BravaisClass::mP),
    ArithmeticCrystalClassEntry::new(8, "2/mC", GeometricCrystalClass::C2h, BravaisClass::mC),
    // Crystal system: orthorhombic
    ArithmeticCrystalClassEntry::new(9, "222P", GeometricCrystalClass::D2, BravaisClass::oP),
    ArithmeticCrystalClassEntry::new(10, "222C", GeometricCrystalClass::D2, BravaisClass::oS),
    ArithmeticCrystalClassEntry::new(11, "222F", GeometricCrystalClass::D2, BravaisClass::oF),
    ArithmeticCrystalClassEntry::new(12, "222I", GeometricCrystalClass::D2, BravaisClass::oI),
    ArithmeticCrystalClassEntry::new(13, "mm2P", GeometricCrystalClass::C2v, BravaisClass::oP),
    ArithmeticCrystalClassEntry::new(14, "mm2C", GeometricCrystalClass::C2v, BravaisClass::oS),
    ArithmeticCrystalClassEntry::new(15, "2mmC", GeometricCrystalClass::C2v, BravaisClass::oS),
    ArithmeticCrystalClassEntry::new(16, "mm2F", GeometricCrystalClass::C2v, BravaisClass::oF),
    ArithmeticCrystalClassEntry::new(17, "mm2I", GeometricCrystalClass::C2v, BravaisClass::oI),
    ArithmeticCrystalClassEntry::new(18, "mmmP", GeometricCrystalClass::D2h, BravaisClass::oP),
    ArithmeticCrystalClassEntry::new(19, "mmmC", GeometricCrystalClass::D2h, BravaisClass::oS),
    ArithmeticCrystalClassEntry::new(20, "mmmF", GeometricCrystalClass::D2h, BravaisClass::oF),
    ArithmeticCrystalClassEntry::new(21, "mmmI", GeometricCrystalClass::D2h, BravaisClass::oI),
    // Crystal system: tetragonal
    ArithmeticCrystalClassEntry::new(22, "4P", GeometricCrystalClass::C4, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(23, "4I", GeometricCrystalClass::C4, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(24, "-4P", GeometricCrystalClass::S4, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(25, "-4I", GeometricCrystalClass::S4, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(26, "4/mP", GeometricCrystalClass::C4h, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(27, "4/mI", GeometricCrystalClass::C4h, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(28, "422P", GeometricCrystalClass::D4, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(29, "422I", GeometricCrystalClass::D4, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(30, "4mmP", GeometricCrystalClass::C4v, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(31, "4mmI", GeometricCrystalClass::C4v, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(32, "-42mP", GeometricCrystalClass::D2d, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(33, "-4m2P", GeometricCrystalClass::D2d, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(34, "-4m2I", GeometricCrystalClass::D2d, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(35, "-42mI", GeometricCrystalClass::D2d, BravaisClass::tI),
    ArithmeticCrystalClassEntry::new(36, "4/mmmP", GeometricCrystalClass::D4h, BravaisClass::tP),
    ArithmeticCrystalClassEntry::new(37, "4/mmmI", GeometricCrystalClass::D4h, BravaisClass::tI),
    // Crystal system: trigonal
    ArithmeticCrystalClassEntry::new(38, "3P", GeometricCrystalClass::C3, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(39, "3R", GeometricCrystalClass::C3, BravaisClass::hR),
    ArithmeticCrystalClassEntry::new(40, "-3P", GeometricCrystalClass::C3i, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(41, "-3R", GeometricCrystalClass::C3i, BravaisClass::hR),
    ArithmeticCrystalClassEntry::new(42, "312P", GeometricCrystalClass::D3, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(43, "321P", GeometricCrystalClass::D3, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(44, "32R", GeometricCrystalClass::D3, BravaisClass::hR),
    ArithmeticCrystalClassEntry::new(45, "3m1P", GeometricCrystalClass::C3v, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(46, "31mP", GeometricCrystalClass::C3v, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(47, "3mR", GeometricCrystalClass::C3v, BravaisClass::hR),
    ArithmeticCrystalClassEntry::new(48, "-31mP", GeometricCrystalClass::D3d, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(49, "-3m1P", GeometricCrystalClass::D3d, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(50, "-3mR", GeometricCrystalClass::D3d, BravaisClass::hR),
    // Crystal system: hexagonal
    ArithmeticCrystalClassEntry::new(51, "6P", GeometricCrystalClass::C6, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(52, "-6P", GeometricCrystalClass::C3h, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(53, "6/mP", GeometricCrystalClass::C6h, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(54, "622P", GeometricCrystalClass::D6, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(55, "6mmP", GeometricCrystalClass::C6v, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(56, "-62mP", GeometricCrystalClass::D3h, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(57, "-6m2P", GeometricCrystalClass::D3h, BravaisClass::hP),
    ArithmeticCrystalClassEntry::new(58, "6/mmmP", GeometricCrystalClass::D6h, BravaisClass::hP),
    // Crystal system: cubic
    ArithmeticCrystalClassEntry::new(59, "23P", GeometricCrystalClass::T, BravaisClass::cP),
    ArithmeticCrystalClassEntry::new(60, "23F", GeometricCrystalClass::T, BravaisClass::cF),
    ArithmeticCrystalClassEntry::new(61, "23I", GeometricCrystalClass::T, BravaisClass::cI),
    ArithmeticCrystalClassEntry::new(62, "m-3P", GeometricCrystalClass::Th, BravaisClass::cP),
    ArithmeticCrystalClassEntry::new(63, "m-3F", GeometricCrystalClass::Th, BravaisClass::cF),
    ArithmeticCrystalClassEntry::new(64, "m-3I", GeometricCrystalClass::Th, BravaisClass::cI),
    ArithmeticCrystalClassEntry::new(65, "432P", GeometricCrystalClass::O, BravaisClass::cP),
    ArithmeticCrystalClassEntry::new(66, "432F", GeometricCrystalClass::O, BravaisClass::cF),
    ArithmeticCrystalClassEntry::new(67, "432I", GeometricCrystalClass::O, BravaisClass::cI),
    ArithmeticCrystalClassEntry::new(68, "-43mP", GeometricCrystalClass::Td, BravaisClass::cP),
    ArithmeticCrystalClassEntry::new(69, "-43mF", GeometricCrystalClass::Td, BravaisClass::cF),
    ArithmeticCrystalClassEntry::new(70, "-43mI", GeometricCrystalClass::Td, BravaisClass::cI),
    ArithmeticCrystalClassEntry::new(71, "m-3mP", GeometricCrystalClass::Oh, BravaisClass::cP),
    ArithmeticCrystalClassEntry::new(72, "m-3mF", GeometricCrystalClass::Oh, BravaisClass::cF),
    ArithmeticCrystalClassEntry::new(73, "m-3mI", GeometricCrystalClass::Oh, BravaisClass::cI),
];
