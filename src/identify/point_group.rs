use std::collections::{HashSet, VecDeque};

use nalgebra::Vector3;

use crate::base::lattice::Lattice;
use crate::base::operation::Rotation;
use crate::base::transformation::TransformationMatrix;
use crate::data::arithmetic_crystal_class::ArithmeticNumber;
use crate::data::classification::{GeometricCrystalClass, LaueClass};
use crate::data::hall_symbol::HallSymbol;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum RotationType {
    Rotation1,      // 1
    Rotation2,      // 2
    Rotation3,      // 3
    Rotation4,      // 4
    Rotation6,      // 6
    RotoInversion1, // -1 = S2
    RotoInversion2, // -2 = m = S1
    RotoInversion3, // -3 = S6^-1
    RotoInversion4, // -4 = S4^-1
    RotoInversion6, // -6 = S3^-1
}

/// Crystallographic point group
#[derive(Debug)]
pub struct PointGroup {
    pub lattice: Lattice,
    pub rotations: Vec<Rotation>,
    //
    pub generators: Vec<usize>,
}

/// Crystallographic point group
impl PointGroup {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>) -> Self {
        let generators = choose_generators(&rotations);
        Self {
            lattice,
            rotations,
            generators,
        }
    }

    /// Construct representative point group from geometric crystal class
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
        let hall_symbol = HallSymbol::from_hall_number(hall_number);
        let operations = hall_symbol.traverse();
        Self::new(operations.lattice, operations.rotations)
    }

    pub fn order(&self) -> usize {
        self.rotations.len()
    }

    /// Let P be the (unimodular) transformation matrix
    /// P^-1 * self.rotations * P = PointGroup::from_arithmetic_crystal_class(arithmetic_number)
    ///
    /// Implement algorithm presented in "R. W. Grosse-Kunstleve. Algorithms for deriving crystallographic space-group information. Acta Cryst. A, 55, 383-395 (1999)".
    pub fn match_arithmetic(&self) -> (ArithmeticNumber, TransformationMatrix) {
        let rotation_types = self
            .rotations
            .iter()
            .map(|rotation| identify_rotation_type(rotation))
            .collect();
        let geometric_crystal_class = identify_geometric_crystal_class(&rotation_types);
        let laue_class = LaueClass::from_geometric_crystal_class(geometric_crystal_class);
        unimplemented!()
    }
}

fn identify_rotation_type(rotation: &Rotation) -> RotationType {
    let tr = rotation.trace();
    let det = rotation.map(|e| e as f64).determinant().round() as i32;

    match (tr, det) {
        (3, 1) => RotationType::Rotation1,
        (-1, 1) => RotationType::Rotation2,
        (0, 1) => RotationType::Rotation3,
        (1, 1) => RotationType::Rotation4,
        (2, 1) => RotationType::Rotation6,
        (-3, -1) => RotationType::RotoInversion1,
        (1, -1) => RotationType::RotoInversion2,
        (0, -1) => RotationType::RotoInversion3,
        (-1, -1) => RotationType::RotoInversion4,
        (-2, -1) => RotationType::RotoInversion6,
        _ => unreachable!("Unknown rotation type"),
    }
}

/// Algorithm in Section 4(a) of GK99
fn find_axis_direction(rotation: &Rotation) -> Vector3<i32> {
    let proper_rotation = if rotation.map(|e| e as f64).determinant().round() < 0.0 {
        rotation.clone()
    } else {
        -1 * rotation.clone()
    };
    let rotation_type = identify_rotation_type(&proper_rotation);

    // Find axis direction by solving eigenvalue problem with eigenvalue 1

    unimplemented!()
}

/// Use look up table in Table 6 of https://arxiv.org/pdf/1808.01590.pdf
fn identify_geometric_crystal_class(rotation_types: &Vec<RotationType>) -> GeometricCrystalClass {
    // count RotationTypes in point_group
    let mut rotation_types_count = [0; 10];
    for rotation_type in rotation_types {
        match rotation_type {
            RotationType::RotoInversion6 => rotation_types_count[0] += 1,
            RotationType::RotoInversion4 => rotation_types_count[1] += 1,
            RotationType::RotoInversion3 => rotation_types_count[2] += 1,
            RotationType::RotoInversion2 => rotation_types_count[3] += 1,
            RotationType::RotoInversion1 => rotation_types_count[4] += 1,
            RotationType::Rotation1 => rotation_types_count[5] += 1,
            RotationType::Rotation2 => rotation_types_count[6] += 1,
            RotationType::Rotation3 => rotation_types_count[7] += 1,
            RotationType::Rotation4 => rotation_types_count[8] += 1,
            RotationType::Rotation6 => rotation_types_count[9] += 1,
        }
    }
    match rotation_types_count {
        // Triclinic
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0] => GeometricCrystalClass::C1,
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0] => GeometricCrystalClass::Ci,
        // Monoclinic
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0] => GeometricCrystalClass::C2,
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0] => GeometricCrystalClass::C1h,
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0] => GeometricCrystalClass::C2h,
        // Orthorhombic
        [0, 0, 0, 0, 0, 1, 3, 0, 0, 0] => GeometricCrystalClass::D2,
        [0, 0, 0, 2, 0, 1, 1, 0, 0, 0] => GeometricCrystalClass::C2v,
        [0, 0, 0, 3, 1, 1, 3, 0, 0, 0] => GeometricCrystalClass::D2h,
        // Tetragonal
        [0, 0, 0, 0, 0, 1, 1, 0, 2, 0] => GeometricCrystalClass::C4,
        [0, 2, 0, 0, 0, 1, 1, 0, 0, 0] => GeometricCrystalClass::S4,
        [0, 2, 0, 1, 1, 1, 1, 0, 2, 0] => GeometricCrystalClass::C4h,
        [0, 0, 0, 0, 0, 1, 5, 0, 2, 0] => GeometricCrystalClass::D4,
        [0, 0, 0, 4, 0, 1, 1, 0, 2, 0] => GeometricCrystalClass::C4v,
        [0, 2, 0, 2, 0, 1, 3, 0, 0, 0] => GeometricCrystalClass::D2d,
        [0, 2, 0, 5, 1, 1, 5, 0, 2, 0] => GeometricCrystalClass::D4h,
        // Trigonal
        [0, 0, 0, 0, 0, 1, 0, 2, 0, 0] => GeometricCrystalClass::C3,
        [0, 0, 2, 0, 1, 1, 0, 2, 0, 0] => GeometricCrystalClass::C3i,
        [0, 0, 0, 0, 0, 1, 3, 2, 0, 0] => GeometricCrystalClass::D3,
        [0, 0, 0, 3, 0, 1, 0, 2, 0, 0] => GeometricCrystalClass::C3v,
        [0, 0, 2, 3, 1, 1, 3, 2, 0, 0] => GeometricCrystalClass::D3d,
        // Hexagonal
        [0, 0, 0, 0, 0, 1, 1, 2, 0, 2] => GeometricCrystalClass::C6,
        [2, 0, 0, 1, 0, 1, 0, 2, 0, 0] => GeometricCrystalClass::C3h,
        [2, 0, 2, 1, 1, 1, 1, 2, 0, 2] => GeometricCrystalClass::C6h,
        [0, 0, 0, 0, 0, 1, 7, 2, 0, 2] => GeometricCrystalClass::D6,
        [0, 0, 0, 6, 0, 1, 1, 2, 0, 2] => GeometricCrystalClass::C6v,
        [2, 0, 0, 4, 0, 1, 3, 2, 0, 0] => GeometricCrystalClass::D3h,
        [2, 0, 2, 7, 1, 1, 7, 2, 0, 2] => GeometricCrystalClass::D6h,
        // Cubic
        [0, 0, 0, 0, 0, 1, 3, 8, 0, 0] => GeometricCrystalClass::T,
        [0, 0, 8, 3, 1, 1, 3, 8, 0, 0] => GeometricCrystalClass::Th,
        [0, 0, 0, 0, 0, 1, 9, 8, 6, 0] => GeometricCrystalClass::O,
        [0, 6, 0, 6, 0, 1, 3, 8, 0, 0] => GeometricCrystalClass::Td,
        [0, 6, 8, 9, 1, 1, 9, 8, 6, 0] => GeometricCrystalClass::Oh,
        // Unknown
        _ => unreachable!("Unknown point group"),
    }
}

fn choose_generators(elements: &Vec<Rotation>) -> Vec<usize> {
    let mut generators = Vec::new();

    let mut visited = HashSet::new();
    let identity = Rotation::identity();
    visited.insert(identity);

    for (i, element) in elements.iter().enumerate() {
        if visited.contains(element) {
            continue;
        }
        generators.push(i);

        let mut queue = VecDeque::new();
        for other in visited.iter().cloned() {
            queue.push_back(other);
        }

        while !queue.is_empty() {
            let other = queue.pop_front().unwrap();
            let product = element * other;
            if !visited.contains(&product) {
                visited.insert(product.clone());
                queue.push_back(product.clone());
            }
        }
    }

    generators
}

#[cfg(test)]
mod tests {
    use std::collections::{HashSet, VecDeque};

    use rstest::rstest;

    use super::PointGroup;
    use crate::base::operation::Rotation;
    use crate::data::classification::GeometricCrystalClass;

    fn traverse(generators: &Vec<Rotation>) -> Vec<Rotation> {
        let mut queue = VecDeque::new();
        let mut group = HashSet::new();

        queue.push_back(Rotation::identity());

        while !queue.is_empty() {
            let element = queue.pop_front().unwrap();
            if group.contains(&element) {
                continue;
            }
            group.insert(element.clone());
            for generator in generators {
                let product = element * generator;
                queue.push_back(product);
            }
        }

        group.into_iter().collect()
    }

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
        let point_group = PointGroup::from_geometric_crystal_class(geometric_crystal_class);
        assert_eq!(point_group.order(), order);

        let elements = traverse(&point_group.rotations);
        assert_eq!(elements.len(), order);
    }
}
