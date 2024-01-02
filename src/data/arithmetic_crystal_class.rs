use crate::data::classification::{BravaisClass, GeometricCrystalClass};

pub type ArithmeticNumber = i32;
pub type ArithmeticCrystalClassSymbol = &'static str;

pub type ArithmeticCrystalClassEntry = (
    ArithmeticNumber,
    ArithmeticCrystalClassSymbol,
    GeometricCrystalClass,
    BravaisClass,
);

// Ordered the same as https://dictionary.iucr.org/Arithmetic_crystal_class
pub const ARITHMETIC_CRYSTAL_CLASS_DATABASE: [ArithmeticCrystalClassEntry; 73] = [
    // Crystal system: triclinic
    (1, "1P", GeometricCrystalClass::C1, BravaisClass::aP),
    (2, "-1P", GeometricCrystalClass::Ci, BravaisClass::aP),
    // Crystal system: monoclinic
    (3, "2P", GeometricCrystalClass::C2, BravaisClass::mP),
    (4, "2C", GeometricCrystalClass::C2, BravaisClass::mC),
    (5, "mP", GeometricCrystalClass::C1h, BravaisClass::mP),
    (6, "mC", GeometricCrystalClass::C1h, BravaisClass::mC),
    (7, "2/mP", GeometricCrystalClass::C2h, BravaisClass::mP),
    (8, "2/mC", GeometricCrystalClass::C2h, BravaisClass::mC),
    // Crystal system: orthorhombic
    (9, "222P", GeometricCrystalClass::D2, BravaisClass::oP),
    (10, "222C", GeometricCrystalClass::D2, BravaisClass::oS),
    (11, "222F", GeometricCrystalClass::D2, BravaisClass::oF),
    (12, "222I", GeometricCrystalClass::D2, BravaisClass::oI),
    (13, "mm2P", GeometricCrystalClass::C2v, BravaisClass::oP),
    (14, "mm2C", GeometricCrystalClass::C2v, BravaisClass::oS),
    (15, "2mmC", GeometricCrystalClass::C2v, BravaisClass::oS),
    (16, "mm2F", GeometricCrystalClass::C2v, BravaisClass::oF),
    (17, "mm2I", GeometricCrystalClass::C2v, BravaisClass::oI),
    (18, "mmmP", GeometricCrystalClass::D2h, BravaisClass::oP),
    (19, "mmmC", GeometricCrystalClass::D2h, BravaisClass::oS),
    (20, "mmmF", GeometricCrystalClass::D2h, BravaisClass::oF),
    (21, "mmmI", GeometricCrystalClass::D2h, BravaisClass::oI),
    // Crystal system: tetragonal
    (22, "4P", GeometricCrystalClass::C4, BravaisClass::tP),
    (23, "4I", GeometricCrystalClass::C4, BravaisClass::tI),
    (24, "-4P", GeometricCrystalClass::S4, BravaisClass::tP),
    (25, "-4I", GeometricCrystalClass::S4, BravaisClass::tI),
    (26, "4/mP", GeometricCrystalClass::C4h, BravaisClass::tP),
    (27, "4/mI", GeometricCrystalClass::C4h, BravaisClass::tI),
    (28, "422P", GeometricCrystalClass::D4, BravaisClass::tP),
    (29, "422I", GeometricCrystalClass::D4, BravaisClass::tI),
    (30, "4mmP", GeometricCrystalClass::C4v, BravaisClass::tP),
    (31, "4mmI", GeometricCrystalClass::C4v, BravaisClass::tI),
    (32, "-42mP", GeometricCrystalClass::D2d, BravaisClass::tP),
    (33, "-4m2P", GeometricCrystalClass::D2d, BravaisClass::tP),
    (34, "-4m2I", GeometricCrystalClass::D2d, BravaisClass::tI),
    (35, "-42mI", GeometricCrystalClass::D2d, BravaisClass::tI),
    (36, "4/mmmP", GeometricCrystalClass::D4h, BravaisClass::tP),
    (37, "4/mmmI", GeometricCrystalClass::D4h, BravaisClass::tI),
    // Crystal system: trigonal
    (38, "3P", GeometricCrystalClass::C3, BravaisClass::hP),
    (39, "3R", GeometricCrystalClass::C3, BravaisClass::hR),
    (40, "-3P", GeometricCrystalClass::C3i, BravaisClass::hP),
    (41, "-3R", GeometricCrystalClass::C3i, BravaisClass::hR),
    (42, "312P", GeometricCrystalClass::D3, BravaisClass::hP),
    (43, "321P", GeometricCrystalClass::D3, BravaisClass::hP),
    (44, "32R", GeometricCrystalClass::D3, BravaisClass::hR),
    (45, "3m1P", GeometricCrystalClass::C3v, BravaisClass::hP),
    (46, "31mP", GeometricCrystalClass::C3v, BravaisClass::hP),
    (47, "3mR", GeometricCrystalClass::C3v, BravaisClass::hR),
    (48, "-31mP", GeometricCrystalClass::D3d, BravaisClass::hP),
    (49, "-3m1P", GeometricCrystalClass::D3d, BravaisClass::hP),
    (50, "-3mR", GeometricCrystalClass::D3d, BravaisClass::hR),
    // Crystal system: hexagonal
    (51, "6P", GeometricCrystalClass::C6, BravaisClass::hP),
    (52, "-6P", GeometricCrystalClass::C3h, BravaisClass::hP),
    (53, "6/mP", GeometricCrystalClass::C6h, BravaisClass::hP),
    (54, "622P", GeometricCrystalClass::D6, BravaisClass::hP),
    (55, "6mmP", GeometricCrystalClass::C6v, BravaisClass::hP),
    (56, "-62mP", GeometricCrystalClass::D3h, BravaisClass::hP),
    (57, "-6m2P", GeometricCrystalClass::D3h, BravaisClass::hP),
    (58, "6/mmmP", GeometricCrystalClass::D6h, BravaisClass::hP),
    // Crystal system: cubic
    (59, "23P", GeometricCrystalClass::T, BravaisClass::cP),
    (60, "23F", GeometricCrystalClass::T, BravaisClass::cF),
    (61, "23I", GeometricCrystalClass::T, BravaisClass::cI),
    (62, "m-3P", GeometricCrystalClass::Th, BravaisClass::cP),
    (63, "m-3F", GeometricCrystalClass::Th, BravaisClass::cF),
    (64, "m-3I", GeometricCrystalClass::Th, BravaisClass::cI),
    (65, "432P", GeometricCrystalClass::O, BravaisClass::cP),
    (66, "432F", GeometricCrystalClass::O, BravaisClass::cF),
    (67, "432I", GeometricCrystalClass::O, BravaisClass::cI),
    (68, "-43mP", GeometricCrystalClass::Td, BravaisClass::cP),
    (69, "-43mF", GeometricCrystalClass::Td, BravaisClass::cF),
    (70, "-43mI", GeometricCrystalClass::Td, BravaisClass::cI),
    (71, "m-3mP", GeometricCrystalClass::Oh, BravaisClass::cP),
    (72, "m-3mF", GeometricCrystalClass::Oh, BravaisClass::cF),
    (73, "m-3mI", GeometricCrystalClass::Oh, BravaisClass::cI),
];

pub fn arithmetic_crystal_class_entry(
    arithmetic_number: ArithmeticNumber,
) -> ArithmeticCrystalClassEntry {
    ARITHMETIC_CRYSTAL_CLASS_DATABASE[arithmetic_number as usize - 1]
}
