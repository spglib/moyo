use crate::base::lattice::Lattice;
use crate::base::operation::{Rotation, RotationType};
use crate::data::classification::GeometricCrystalClass;
use crate::data::hall_symbol::HallSymbol;

/// Crystallographic point group
#[derive(Debug)]
pub struct PointGroup {
    pub lattice: Lattice,
    pub rotations: Vec<Rotation>,
}

impl PointGroup {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>) -> Self {
        Self { lattice, rotations }
    }

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
        let hall_symbol = HallSymbol::from_hall_number(hall_number);
        let operations = hall_symbol.traverse();
        Self::new(operations.lattice, operations.rotations)
    }

    pub fn order(&self) -> usize {
        self.rotations.len()
    }
}

fn identify_rotation_type(rotation: &Rotation) -> Option<RotationType> {
    let tr = rotation.trace();
    let det = rotation.map(|e| e as f64).determinant().round() as i32;

    match (tr, det) {
        (3, 1) => Some(RotationType::Rotation1),
        (-1, 1) => Some(RotationType::Rotation2),
        (0, 1) => Some(RotationType::Rotation3),
        (1, 1) => Some(RotationType::Rotation4),
        (2, 1) => Some(RotationType::Rotation6),
        (-3, -1) => Some(RotationType::RotoInversion1),
        (1, -1) => Some(RotationType::RotoInversion2),
        (0, -1) => Some(RotationType::RotoInversion3),
        (-1, -1) => Some(RotationType::RotoInversion4),
        (-2, -1) => Some(RotationType::RotoInversion6),
        _ => None,
    }
}

/// Use look up table in Table 6 of https://arxiv.org/pdf/1808.01590.pdf
pub fn identify_geometric_crystal_class(point_group: &PointGroup) -> Option<GeometricCrystalClass> {
    // count RotationTypes in point_group
    let mut rotation_types = [0; 10];
    for rotation in &point_group.rotations {
        if let Some(rotation_type) = identify_rotation_type(rotation) {
            match rotation_type {
                RotationType::RotoInversion6 => rotation_types[0] += 1,
                RotationType::RotoInversion4 => rotation_types[1] += 1,
                RotationType::RotoInversion3 => rotation_types[2] += 1,
                RotationType::RotoInversion2 => rotation_types[3] += 1,
                RotationType::RotoInversion1 => rotation_types[4] += 1,
                RotationType::Rotation1 => rotation_types[5] += 1,
                RotationType::Rotation2 => rotation_types[6] += 1,
                RotationType::Rotation3 => rotation_types[7] += 1,
                RotationType::Rotation4 => rotation_types[8] += 1,
                RotationType::Rotation6 => rotation_types[9] += 1,
            }
        } else {
            return None;
        }
    }
    match rotation_types {
        // Triclinic
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::C1),
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::Ci),
        // Monoclinic
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2),
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::C1h),
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2h),
        // Orthorhombic
        [0, 0, 0, 0, 0, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2),
        [0, 0, 0, 2, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2v),
        [0, 0, 0, 3, 1, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2h),
        // Tetragonal
        [0, 0, 0, 0, 0, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4),
        [0, 2, 0, 0, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::S4),
        [0, 2, 0, 1, 1, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4h),
        [0, 0, 0, 0, 0, 1, 5, 0, 2, 0] => Some(GeometricCrystalClass::D4),
        [0, 0, 0, 4, 0, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4v),
        [0, 2, 0, 2, 0, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2d),
        [0, 2, 0, 5, 1, 1, 5, 0, 2, 0] => Some(GeometricCrystalClass::D4h),
        // Trigonal
        [0, 0, 0, 0, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3),
        [0, 0, 2, 0, 1, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3i),
        [0, 0, 0, 0, 0, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3),
        [0, 0, 0, 3, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3v),
        [0, 0, 2, 3, 1, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3d),
        // Hexagonal
        [0, 0, 0, 0, 0, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6),
        [2, 0, 0, 1, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3h),
        [2, 0, 2, 1, 1, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6h),
        [0, 0, 0, 0, 0, 1, 7, 2, 0, 2] => Some(GeometricCrystalClass::D6),
        [0, 0, 0, 6, 0, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6v),
        [2, 0, 0, 4, 0, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3h),
        [2, 0, 2, 7, 1, 1, 7, 2, 0, 2] => Some(GeometricCrystalClass::D6h),
        // Cubic
        [0, 0, 0, 0, 0, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::T),
        [0, 0, 8, 3, 1, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::Th),
        [0, 0, 0, 0, 0, 1, 9, 8, 6, 0] => Some(GeometricCrystalClass::O),
        [0, 6, 0, 6, 0, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::Td),
        [0, 6, 8, 9, 1, 1, 9, 8, 6, 0] => Some(GeometricCrystalClass::Oh),
        // Unknown
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::PointGroup;
    use crate::data::classification::GeometricCrystalClass;

    #[rstest]
    // Triclinic
    #[case(GeometricCrystalClass::C1, 1)]
    #[case(GeometricCrystalClass::Ci, 2)]
    // Monoclinic
    #[case(GeometricCrystalClass::C2, 2)]
    #[case(GeometricCrystalClass::C1h, 2)]
    #[case(GeometricCrystalClass::C2h, 4)]
    // Orthorhombic
    #[case(GeometricCrystalClass::D2, 4)]
    #[case(GeometricCrystalClass::C2v, 4)]
    #[case(GeometricCrystalClass::D2h, 8)]
    // Tetragonal
    #[case(GeometricCrystalClass::C4, 4)]
    #[case(GeometricCrystalClass::S4, 4)]
    #[case(GeometricCrystalClass::C4h, 8)]
    #[case(GeometricCrystalClass::D4, 8)]
    #[case(GeometricCrystalClass::C4v, 8)]
    #[case(GeometricCrystalClass::D2d, 8)]
    #[case(GeometricCrystalClass::D4h, 16)]
    // Trigonal
    #[case(GeometricCrystalClass::C3, 3)]
    #[case(GeometricCrystalClass::C3i, 6)]
    #[case(GeometricCrystalClass::D3, 6)]
    #[case(GeometricCrystalClass::C3v, 6)]
    #[case(GeometricCrystalClass::D3d, 12)]
    // Hexagonal
    #[case(GeometricCrystalClass::C6, 6)]
    #[case(GeometricCrystalClass::C3h, 6)]
    #[case(GeometricCrystalClass::C6h, 12)]
    #[case(GeometricCrystalClass::D6, 12)]
    #[case(GeometricCrystalClass::C6v, 12)]
    #[case(GeometricCrystalClass::D3h, 12)]
    #[case(GeometricCrystalClass::D6h, 24)]
    // Cubic
    #[case(GeometricCrystalClass::T, 12)]
    #[case(GeometricCrystalClass::Th, 24)]
    #[case(GeometricCrystalClass::O, 24)]
    #[case(GeometricCrystalClass::Td, 24)]
    #[case(GeometricCrystalClass::Oh, 48)]
    fn test_point_group_representative(
        #[case] geometric_crystal_class: GeometricCrystalClass,
        #[case] order: usize,
    ) {
        dbg!(&geometric_crystal_class);
        let point_group = PointGroup::from_geometric_crystal_class(geometric_crystal_class);
        assert_eq!(point_group.order(), order);
    }
}
