use std::collections::{HashSet, VecDeque};

use itertools::Itertools;
use nalgebra::{Dyn, OMatrix, OVector, U9};

use crate::base::error::MoyoError;
use crate::base::operation::Rotation;
use crate::base::transformation::TransformationMatrix;
use crate::data::arithmetic_crystal_class::{ArithmeticNumber, ARITHMETIC_CRYSTAL_CLASS_DATABASE};
use crate::data::classification::{CrystalSystem, GeometricCrystalClass};
use crate::data::hall_symbol::{Centering, HallSymbol};
use crate::math::integer_system::IntegerLinearSystem;

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

impl RotationType {
    pub fn value(&self) -> i32 {
        match self {
            RotationType::Rotation1 => 1,
            RotationType::Rotation2 => 2,
            RotationType::Rotation3 => 3,
            RotationType::Rotation4 => 4,
            RotationType::Rotation6 => 6,
            RotationType::RotoInversion1 => -1,
            RotationType::RotoInversion2 => -2,
            RotationType::RotoInversion3 => -3,
            RotationType::RotoInversion4 => -4,
            RotationType::RotoInversion6 => -6,
        }
    }
}

#[derive(Debug)]
/// Specific crystallographic point group in database
pub struct PointGroupRepresentative {
    pub generators: Vec<Rotation>,
    pub centering: Centering,
}

impl PointGroupRepresentative {
    fn new(generators: Vec<Rotation>, centering: Centering) -> Self {
        Self {
            generators,
            centering,
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
        Self::new(hall_symbol.generators.rotations, hall_symbol.lattice_symbol)
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
        Self::new(hall_symbol.generators.rotations, hall_symbol.lattice_symbol)
    }

    pub fn primitive_generators(&self) -> Vec<Rotation> {
        let prim_trans_mat_inv = self.centering.transformation_matrix().map(|e| e as f64);
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

/// Crystallographic point group with group-type information
/// Assume the rotations are given in the (Minkowski-reduced) primitive basis
#[derive(Debug)]
pub struct PointGroup {
    pub arithmetic_number: ArithmeticNumber,
    /// Transformation matrix to the representative for `arithmetic_number`
    /// trans_mat^-1 * self.rotations * trans_mat = AbstractPointGroup::from_arithmetic_crystal_class(arithmetic_number)
    pub trans_mat: TransformationMatrix,
}

pub fn identify_point_group(prim_rotations: &Vec<Rotation>) -> Result<PointGroup, MoyoError> {
    let rotation_types = prim_rotations
        .iter()
        .map(|rotation| identify_rotation_type(rotation))
        .collect();
    let geometric_crystal_class = identify_geometric_crystal_class(&rotation_types);

    let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);
    match crystal_system {
        // Skip triclinic cases as trivial
        CrystalSystem::Triclinic => match geometric_crystal_class {
            GeometricCrystalClass::C1 => {
                return Ok(PointGroup {
                    arithmetic_number: 1,
                    trans_mat: TransformationMatrix::identity(),
                })
            }
            GeometricCrystalClass::Ci => {
                return Ok(PointGroup {
                    arithmetic_number: 2,
                    trans_mat: TransformationMatrix::identity(),
                });
            }
            _ => unreachable!(),
        },
        CrystalSystem::Cubic => {
            if let Some((arithmetic_number, trans_mat)) = match_with_cubic_point_group(
                &prim_rotations,
                &rotation_types,
                geometric_crystal_class,
            ) {
                return Ok(PointGroup {
                    arithmetic_number,
                    trans_mat,
                });
            }
        }
        _ => {
            if let Some((arithmetic_number, trans_mat)) =
                match_with_point_group(&prim_rotations, &rotation_types, geometric_crystal_class)
            {
                return Ok(PointGroup {
                    arithmetic_number,
                    trans_mat,
                });
            }
        }
    }
    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

/// Faster matching algorithm for cubic point groups
fn match_with_cubic_point_group(
    prim_rotations: &Vec<Rotation>,
    rotation_types: &Vec<RotationType>,
    geometric_crystal_class: GeometricCrystalClass,
) -> Option<(ArithmeticNumber, TransformationMatrix)> {
    let generators = choose_generators(&prim_rotations);

    let arithmetic_crystal_class_candidates = ARITHMETIC_CRYSTAL_CLASS_DATABASE
        .iter()
        .filter_map(|(arithmetic_number_db, _, geometric_crystal_class_db, _)| {
            if *geometric_crystal_class_db != geometric_crystal_class {
                return None;
            }

            // Check if prim_trans_mat^-1 * rotation * prim_trans_mat is integer matrix (for all rotation in generators)
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(*arithmetic_number_db);
            let centering = point_group_db.centering;
            let prim_trans_mat = centering.transformation_matrix(); // primitive -> conventional
            let prim_trans_mat_inv = centering.inverse();
            let compatible = generators.iter().all(|i| {
                let rotation = prim_rotations[*i];
                let tmp_mat = prim_trans_mat_inv * (rotation * prim_trans_mat).map(|e| e as f64);
                let tmp_mat_int = tmp_mat.map(|e| e.round() as i32);
                (rotation * prim_trans_mat - prim_trans_mat * tmp_mat_int)
                    .iter()
                    .all(|e| *e == 0)
            });
            if !compatible {
                return None;
            }

            Some((*arithmetic_number_db, point_group_db, centering))
        })
        .collect::<Vec<_>>();
    let primitive_arithmetic_crystal_class = arithmetic_crystal_class_candidates
        .iter()
        .filter(|(_, _, centering)| *centering == Centering::P)
        .next()
        .unwrap();

    let order = prim_rotations.len();
    let other_prim_generators = primitive_arithmetic_crystal_class.1.primitive_generators();

    // Try to map generators
    let other_generator_rotation_types = other_prim_generators
        .iter()
        .map(|rotation| identify_rotation_type(rotation))
        .collect::<Vec<_>>();
    let candidates: Vec<Vec<usize>> = other_generator_rotation_types
        .iter()
        .map(|&rotation_type| {
            (0..order)
                .filter(|&i| rotation_types[i] == rotation_type)
                .collect()
        })
        .collect();

    for pivot in candidates
        .iter()
        .map(|v| v.iter())
        .multi_cartesian_product()
    {
        // Solve P^-1 * self.rotations[self.rotations[pivot[i]]] * P = other.generators[i] (for all i)
        // vec(A * P - P * B) = (I_3 \otimes A - B^T \otimes I_3) * vec(P)
        let mut coeffs = OMatrix::<i32, Dyn, U9>::zeros(9 * other_prim_generators.len());
        for (k, &pk) in pivot.iter().enumerate() {
            let rotation = prim_rotations[*pk];
            let other_rotation = other_prim_generators[k];
            let identity = Rotation::identity();
            let adj =
                identity.kronecker(&rotation) - other_rotation.transpose().kronecker(&identity);
            for i in 0..9 {
                for j in 0..9 {
                    coeffs[(9 * k + i, j)] = adj[(i, j)];
                }
            }
        }
        let solution =
            IntegerLinearSystem::new(&coeffs, &OVector::<i32, Dyn>::zeros(coeffs.nrows()));

        // Search integer linear combination such that the transformation matrix is unimodular
        if let Some(solution) = solution {
            let trans_mat_basis: Vec<_> = solution
                .nullspace
                .row_iter()
                .map(|e| {
                    TransformationMatrix::new(
                        e[0], e[3], e[6], //
                        e[1], e[4], e[7], //
                        e[2], e[5], e[8], //
                    )
                })
                .collect();

            // trans_mat: self -> primitive
            // The dimension of linear integer system should be one for cubic.
            assert_eq!(trans_mat_basis.len(), 1);
            let mut trans_mat = trans_mat_basis[0];

            // Guarantee det > 0
            let mut det = trans_mat.map(|e| e as f64).determinant().round() as i32;
            if det < 0 {
                trans_mat *= -1;
                det *= -1;
            } else if det == 0 {
                continue;
            }

            for (arithmetic_crystal_class, _, centering) in
                arithmetic_crystal_class_candidates.iter()
            {
                if centering.order() as i32 != det {
                    continue;
                }

                // trans_mat_conv: self -> conventional
                let trans_mat_conv = centering.transformation_matrix() * trans_mat;
                return Some((*arithmetic_crystal_class, trans_mat_conv));
            }
        }
    }

    None
}

fn match_with_point_group(
    prim_rotations: &Vec<Rotation>,
    rotation_types: &Vec<RotationType>,
    geometric_crystal_class: GeometricCrystalClass,
) -> Option<(ArithmeticNumber, TransformationMatrix)> {
    let order = prim_rotations.len();

    for (arithmetic_number_db, _, geometric_crystal_class_db, _) in
        ARITHMETIC_CRYSTAL_CLASS_DATABASE.iter()
    {
        if *geometric_crystal_class_db != geometric_crystal_class {
            continue;
        }

        let point_group_db =
            PointGroupRepresentative::from_arithmetic_crystal_class(*arithmetic_number_db);
        let other_prim_generators = point_group_db.primitive_generators();

        // Try to map generators
        let other_generator_rotation_types = other_prim_generators
            .iter()
            .map(|rotation| identify_rotation_type(rotation))
            .collect::<Vec<_>>();
        let candidates: Vec<Vec<usize>> = other_generator_rotation_types
            .iter()
            .map(|&rotation_type| {
                (0..order)
                    .filter(|&i| rotation_types[i] == rotation_type)
                    .collect()
            })
            .collect();

        for pivot in candidates
            .iter()
            .map(|v| v.iter())
            .multi_cartesian_product()
        {
            // Solve self.rotations[self.rotations[pivot[i]]] * P = P * other.generators[i] (for all i)
            // vec(A * P - P * B) = (I_3 \otimes A - B^T \otimes I_3) * vec(P)
            let mut coeffs = OMatrix::<i32, Dyn, U9>::zeros(9 * other_prim_generators.len());
            for (k, &pk) in pivot.iter().enumerate() {
                let rotation = prim_rotations[*pk];
                let other_rotation = other_prim_generators[k];
                let identity = Rotation::identity();
                let adj =
                    identity.kronecker(&rotation) - other_rotation.transpose().kronecker(&identity);
                for i in 0..9 {
                    for j in 0..9 {
                        coeffs[(9 * k + i, j)] = adj[(i, j)];
                    }
                }
            }
            let solution =
                IntegerLinearSystem::new(&coeffs, &OVector::<i32, Dyn>::zeros(coeffs.nrows()));

            // Search integer linear combination such that the transformation matrix is unimodular
            if let Some(solution) = solution {
                let trans_mat_basis: Vec<_> = solution
                    .nullspace
                    .row_iter()
                    .map(|e| {
                        TransformationMatrix::new(
                            e[0], e[3], e[6], //
                            e[1], e[4], e[7], //
                            e[2], e[5], e[8], //
                        )
                    })
                    .collect();

                // Consider coefficients in [-2, 2], which will be sufficient for Delaunay reduced basis
                for comb in (0..trans_mat_basis.len())
                    .map(|_| -2..=2)
                    .multi_cartesian_product()
                {
                    let mut trans_mat = TransformationMatrix::zeros();
                    for (i, matrix) in trans_mat_basis.iter().enumerate() {
                        trans_mat += comb[i] * matrix;
                    }
                    let det = trans_mat.map(|e| e as f64).determinant().round() as i32;
                    if det == 1 {
                        return Some((*arithmetic_number_db, trans_mat));
                    }
                }
            }
        }
    }

    None
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

    use strum::IntoEnumIterator;

    use super::{identify_point_group, PointGroupRepresentative};
    use crate::base::operation::Rotation;
    use crate::data::classification::GeometricCrystalClass;

    fn traverse(generators: &Vec<Rotation>) -> Vec<Rotation> {
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();
        let mut group = vec![];

        queue.push_back(Rotation::identity());

        while !queue.is_empty() {
            let element = queue.pop_front().unwrap();
            if visited.contains(&element) {
                continue;
            }
            visited.insert(element.clone());
            group.push(element.clone());

            for generator in generators {
                let product = element * generator;
                queue.push_back(product);
            }
        }

        group
    }

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

    #[test]
    fn test_point_group_match() {
        // TODO: this test takes ~6 seconds with debug build. We may need to speed up `PointGroup::match_arithmetic`
        for arithmetic_number in 1..=73 {
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
            let primitive_generators = point_group_db.primitive_generators();
            let prim_rotations = traverse(&primitive_generators);
            let point_group = identify_point_group(&prim_rotations).unwrap();
            assert_eq!(point_group.arithmetic_number, arithmetic_number);
        }
    }
}
