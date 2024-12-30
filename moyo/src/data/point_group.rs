use super::arithmetic_crystal_class::ArithmeticNumber;
use super::centering::Centering;
use super::classification::GeometricCrystalClass;
use super::hall_symbol::HallSymbol;
use crate::base::{project_rotations, Rotations};

#[derive(Debug)]
/// Specific crystallographic point group in database
pub struct PointGroupRepresentative {
    pub generators: Rotations,
    pub centering: Centering,
}

impl PointGroupRepresentative {
    fn new(generators: Rotations, centering: Centering) -> Self {
        Self {
            generators,
            centering,
        }
    }

    /// Construct representative point group from geometric crystal class
    #[allow(dead_code)]
    pub fn from_geometric_crystal_class(geometric_crystal_class: GeometricCrystalClass) -> Self {
        let hall_number = match geometric_crystal_class {
            // Triclinic
            GeometricCrystalClass::C1 => 1,
            GeometricCrystalClass::Ci => 2,
            // Monoclinic (unique axis b)
            GeometricCrystalClass::C2 => 4,
            GeometricCrystalClass::C1h => 18,
            GeometricCrystalClass::C2h => 57,
            // Orthorhombic
            GeometricCrystalClass::D2 => 108,
            GeometricCrystalClass::C2v => 125, // mm2
            GeometricCrystalClass::D2h => 227,
            // Tetragonal
            GeometricCrystalClass::C4 => 349,
            GeometricCrystalClass::S4 => 355,
            GeometricCrystalClass::C4h => 357,
            GeometricCrystalClass::D4 => 366,
            GeometricCrystalClass::C4v => 376,
            GeometricCrystalClass::D2d => 388, // -42m
            GeometricCrystalClass::D4h => 400,
            // Trigonal
            GeometricCrystalClass::C3 => 430,
            GeometricCrystalClass::C3i => 435,
            GeometricCrystalClass::D3 => 438,  // 312
            GeometricCrystalClass::C3v => 446, // 3m1
            GeometricCrystalClass::D3d => 454, // -31m
            // Hexagonal
            GeometricCrystalClass::C6 => 462,
            GeometricCrystalClass::C3h => 468,
            GeometricCrystalClass::C6h => 469,
            GeometricCrystalClass::D6 => 471,
            GeometricCrystalClass::C6v => 477,
            GeometricCrystalClass::D3h => 481, // -6m2
            GeometricCrystalClass::D6h => 485,
            // Cubic
            GeometricCrystalClass::T => 489,
            GeometricCrystalClass::Th => 494,
            GeometricCrystalClass::O => 503,
            GeometricCrystalClass::Td => 511,
            GeometricCrystalClass::Oh => 517,
        };
        let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        Self::new(
            project_rotations(&hall_symbol.generators),
            hall_symbol.centering,
        )
    }

    pub fn from_arithmetic_crystal_class(arithmetic_number: ArithmeticNumber) -> Self {
        // Choose hexagonal axes for rhombohedral space groups
        let hall_number = match arithmetic_number {
            // Triclinic
            1 => 1,
            2 => 2,
            // Monoclinic (unique axis b, cell choice 1)
            3 => 3,
            4 => 9,
            5 => 18,
            6 => 30,
            7 => 57,
            8 => 63,
            // Orthorhombic (setting abc)
            9 => 108,
            10 => 119, // 222C
            11 => 122,
            12 => 123,
            13 => 125,
            14 => 173,
            15 => 185, // mm2A
            16 => 209,
            17 => 215,
            18 => 227,
            19 => 310, // mmmC
            20 => 334,
            21 => 337,
            // Tetragonal
            22 => 349,
            23 => 353,
            24 => 355,
            25 => 356,
            26 => 357,
            27 => 363,
            28 => 366,
            29 => 374,
            30 => 376,
            31 => 384,
            32 => 388,
            33 => 392,
            34 => 396,
            35 => 398,
            36 => 400,
            37 => 424,
            // Trigonal
            38 => 430,
            39 => 433,
            40 => 435,
            41 => 436,
            42 => 438,
            43 => 439,
            44 => 444,
            45 => 446,
            46 => 447,
            47 => 450,
            48 => 454,
            49 => 456,
            50 => 458,
            // Hexagonal
            51 => 462,
            52 => 468,
            53 => 469,
            54 => 471,
            55 => 477,
            56 => 483,
            57 => 481,
            58 => 485,
            // Cubic
            59 => 489,
            60 => 490,
            61 => 491,
            62 => 494,
            63 => 497,
            64 => 500,
            65 => 503,
            66 => 505,
            67 => 507,
            68 => 511,
            69 => 512,
            70 => 513,
            71 => 517,
            72 => 523,
            73 => 529,
            _ => panic!("Invalid arithmetic number"),
        };
        let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        Self::new(
            project_rotations(&hall_symbol.generators),
            hall_symbol.centering,
        )
    }

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
    use strum::IntoEnumIterator;

    use super::PointGroupRepresentative;
    use crate::base::traverse;
    use crate::data::classification::GeometricCrystalClass;

    fn order(geometric_crystal_class: GeometricCrystalClass) -> usize {
        match geometric_crystal_class {
            // Triclinic
            GeometricCrystalClass::C1 => 1,
            GeometricCrystalClass::Ci => 2,
            // Monoclinic
            GeometricCrystalClass::C2 => 2,
            GeometricCrystalClass::C1h => 2,
            GeometricCrystalClass::C2h => 4,
            // Orthorhombic
            GeometricCrystalClass::D2 => 4,
            GeometricCrystalClass::C2v => 4,
            GeometricCrystalClass::D2h => 8,
            // Tetragonal
            GeometricCrystalClass::C4 => 4,
            GeometricCrystalClass::S4 => 4,
            GeometricCrystalClass::C4h => 8,
            GeometricCrystalClass::D4 => 8,
            GeometricCrystalClass::C4v => 8,
            GeometricCrystalClass::D2d => 8,
            GeometricCrystalClass::D4h => 16,
            // Trigonal
            GeometricCrystalClass::C3 => 3,
            GeometricCrystalClass::C3i => 6,
            GeometricCrystalClass::D3 => 6,
            GeometricCrystalClass::C3v => 6,
            GeometricCrystalClass::D3d => 12,
            // Hexagonal
            GeometricCrystalClass::C6 => 6,
            GeometricCrystalClass::C3h => 6,
            GeometricCrystalClass::C6h => 12,
            GeometricCrystalClass::D6 => 12,
            GeometricCrystalClass::C6v => 12,
            GeometricCrystalClass::D3h => 12,
            GeometricCrystalClass::D6h => 24,
            // Cubic
            GeometricCrystalClass::T => 12,
            GeometricCrystalClass::Td => 24,
            GeometricCrystalClass::O => 24,
            GeometricCrystalClass::Th => 24,
            GeometricCrystalClass::Oh => 48,
        }
    }

    #[test]
    // Triclinic
    fn test_point_group_representative() {
        for geometric_crystal_class in GeometricCrystalClass::iter() {
            let point_group =
                PointGroupRepresentative::from_geometric_crystal_class(geometric_crystal_class);
            let rotations = traverse(&point_group.generators);
            assert_eq!(rotations.len(), order(geometric_crystal_class));
        }
    }
}
