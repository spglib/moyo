use std::cmp::Ordering;
use std::collections::{HashSet, VecDeque};

use itertools::Itertools;
use nalgebra::{Dyn, Matrix3, OMatrix, OVector, U9};

use crate::base::error::MoyoError;
use crate::base::operation::Rotation;
use crate::base::transformation::TransformationMatrix;
use crate::data::arithmetic_crystal_class::{iter_arithmetic_crystal_entry, ArithmeticNumber};
use crate::data::classification::{CrystalSystem, GeometricCrystalClass};
use crate::data::hall_symbol::Centering;
use crate::data::point_group::PointGroupRepresentative;
use crate::math::integer_system::IntegerLinearSystem;

/// Crystallographic point group with group-type information
#[derive(Debug)]
pub struct PointGroup {
    pub arithmetic_number: ArithmeticNumber,
    /// Transformation matrix to the representative for `arithmetic_number` in the primitive basis
    pub prim_trans_mat: TransformationMatrix,
    /// Transformation matrix to the representative for `arithmetic_number` in the conventional basis
    pub conv_trans_mat: TransformationMatrix,
}

impl PointGroup {
    /// Identify the arithmetic crystal class from the given rotations and transformation matrix into the representative
    /// Assume the rotations are given in the (reduced) primitive basis
    pub fn new(prim_rotations: &Vec<Rotation>) -> Result<Self, MoyoError> {
        let rotation_types = prim_rotations.iter().map(identify_rotation_type).collect();
        let geometric_crystal_class = identify_geometric_crystal_class(&rotation_types);

        let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);
        match crystal_system {
            // Skip triclinic cases as trivial
            CrystalSystem::Triclinic => match geometric_crystal_class {
                GeometricCrystalClass::C1 => Ok(PointGroup {
                    arithmetic_number: 1,
                    prim_trans_mat: TransformationMatrix::identity(),
                    conv_trans_mat: TransformationMatrix::identity(),
                }),
                GeometricCrystalClass::Ci => Ok(PointGroup {
                    arithmetic_number: 2,
                    prim_trans_mat: TransformationMatrix::identity(),
                    conv_trans_mat: TransformationMatrix::identity(),
                }),
                _ => unreachable!(),
            },
            CrystalSystem::Cubic => match_with_cubic_point_group(
                prim_rotations,
                &rotation_types,
                geometric_crystal_class,
            ),
            _ => match_with_point_group(prim_rotations, &rotation_types, geometric_crystal_class),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
enum RotationType {
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

/// Faster matching algorithm for cubic point groups
fn match_with_cubic_point_group(
    prim_rotations: &Vec<Rotation>,
    rotation_types: &[RotationType],
    geometric_crystal_class: GeometricCrystalClass,
) -> Result<PointGroup, MoyoError> {
    let generators = choose_generators(prim_rotations);

    let arithmetic_crystal_class_candidates = iter_arithmetic_crystal_entry()
        .filter_map(|(arithmetic_number_db, _, geometric_crystal_class_db, _)| {
            if geometric_crystal_class_db != geometric_crystal_class {
                return None;
            }

            // Check if prim_trans_mat^-1 * rotation * prim_trans_mat is integer matrix (for all rotation in generators)
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number_db);
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

            Some((arithmetic_number_db, point_group_db, centering))
        })
        .collect::<Vec<_>>();
    let primitive_arithmetic_crystal_class = arithmetic_crystal_class_candidates
        .iter()
        .find(|(_, _, centering)| *centering == Centering::P)
        .unwrap();

    let order = prim_rotations.len();
    let other_prim_generators = primitive_arithmetic_crystal_class.1.primitive_generators();

    // Try to map generators
    let other_generator_rotation_types = other_prim_generators
        .iter()
        .map(identify_rotation_type)
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
        if let Some(trans_mat_basis) = sylvester(
            &pivot
                .iter()
                .map(|&i| prim_rotations[*i])
                .collect::<Vec<_>>(),
            &other_prim_generators,
        ) {
            // conv_trans_mat: self -> conventional
            // The dimension of linear integer system should be one for cubic.
            assert_eq!(trans_mat_basis.len(), 1);
            let mut conv_trans_mat = trans_mat_basis[0];

            // Guarantee det > 0
            let mut det = conv_trans_mat.map(|e| e as f64).determinant().round() as i32;
            match det.cmp(&0) {
                Ordering::Less => {
                    conv_trans_mat *= -1;
                    det *= -1;
                }
                Ordering::Equal => continue,
                Ordering::Greater => {}
            }

            for (arithmetic_crystal_class, _, centering) in
                arithmetic_crystal_class_candidates.iter()
            {
                if centering.order() as i32 != det {
                    continue;
                }
                // conv_trans_mat: self -> conventional
                let prim_trans_mat = (conv_trans_mat.map(|e| e as f64) * centering.inverse())
                    .map(|e| e.round() as i32);

                return Ok(PointGroup {
                    arithmetic_number: *arithmetic_crystal_class,
                    prim_trans_mat,
                    conv_trans_mat,
                });
            }
        }
    }

    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

fn match_with_point_group(
    prim_rotations: &Vec<Rotation>,
    rotation_types: &[RotationType],
    geometric_crystal_class: GeometricCrystalClass,
) -> Result<PointGroup, MoyoError> {
    let order = prim_rotations.len();

    for (arithmetic_number_db, _, geometric_crystal_class_db, _) in iter_arithmetic_crystal_entry()
    {
        if geometric_crystal_class_db != geometric_crystal_class {
            continue;
        }

        let point_group_db =
            PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number_db);
        let other_prim_generators = point_group_db.primitive_generators();

        // Try to map generators
        let other_generator_rotation_types = other_prim_generators
            .iter()
            .map(identify_rotation_type)
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
            if let Some(trans_mat_basis) = sylvester(
                &pivot
                    .iter()
                    .map(|&i| prim_rotations[*i])
                    .collect::<Vec<_>>(),
                &other_prim_generators,
            ) {
                // Search integer linear combination such that the transformation matrix is unimodular
                // Consider coefficients in [-2, 2], which will be sufficient for Delaunay reduced basis
                for comb in (0..trans_mat_basis.len())
                    .map(|_| -2..=2)
                    .multi_cartesian_product()
                {
                    // prim_trans_mat: self -> DB(primitive)
                    let mut prim_trans_mat = TransformationMatrix::zeros();
                    for (i, matrix) in trans_mat_basis.iter().enumerate() {
                        prim_trans_mat += comb[i] * matrix;
                    }
                    let det = prim_trans_mat.map(|e| e as f64).determinant().round() as i32;
                    if det == 1 {
                        // conv_trans_mat: self -> DB(conventional)
                        let conv_trans_mat =
                            prim_trans_mat * point_group_db.centering.transformation_matrix();
                        return Ok(PointGroup {
                            arithmetic_number: arithmetic_number_db,
                            prim_trans_mat,
                            conv_trans_mat,
                        });
                    }
                }
            }
        }
    }

    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

/// Solve P^-1 * A[i] * P = B[i] (for all i)
/// vec(A * P - P * B) = (I_3 \otimes A - B^T \otimes I_3) * vec(P)
fn sylvester(a: &Vec<Matrix3<i32>>, b: &Vec<Matrix3<i32>>) -> Option<Vec<Matrix3<i32>>> {
    let size = a.len();
    assert_eq!(size, b.len());

    let mut coeffs = OMatrix::<i32, Dyn, U9>::zeros(9 * size);
    let identity = Rotation::identity();
    for k in 0..size {
        let adj = identity.kronecker(&a[k]) - b[k].transpose().kronecker(&identity);
        for i in 0..9 {
            for j in 0..9 {
                coeffs[(9 * k + i, j)] = adj[(i, j)];
            }
        }
    }
    let solution = IntegerLinearSystem::new(&coeffs, &OVector::<i32, Dyn>::zeros(coeffs.nrows()));

    if let Some(solution) = solution {
        let basis: Vec<_> = solution
            .nullspace
            .row_iter()
            .map(|e| {
                // Vectorization operator is column-major
                TransformationMatrix::new(
                    e[0], e[1], e[2], //
                    e[3], e[4], e[5], //
                    e[6], e[7], e[8], //
                )
                .transpose()
            })
            .collect();
        Some(basis)
    } else {
        None
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

fn choose_generators(elements: &[Rotation]) -> Vec<usize> {
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
                visited.insert(product);
                queue.push_back(product);
            }
        }
    }

    generators
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::{PointGroup, PointGroupRepresentative};
    use crate::base::operation::traverse;

    #[test]
    fn test_point_group_match() {
        // TODO: this test takes ~6 seconds with debug build. We may need to speed up `PointGroup::match_arithmetic`
        for arithmetic_number in 1..=73 {
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
            let primitive_generators = point_group_db.primitive_generators();
            let prim_rotations = traverse(&primitive_generators);
            let point_group = PointGroup::new(&prim_rotations).unwrap();

            assert_eq!(point_group.arithmetic_number, arithmetic_number);
            assert_eq!(
                point_group
                    .prim_trans_mat
                    .map(|e| e as f64)
                    .determinant()
                    .round() as i32,
                1
            );

            // Compare rotations with the primitive
            let mut prim_rotations_set = HashSet::new();
            for rotation in prim_rotations.iter() {
                prim_rotations_set.insert(rotation.clone());
            }
            let prim_trans_inv = point_group
                .prim_trans_mat
                .map(|e| e as f64)
                .try_inverse()
                .unwrap()
                .map(|e| e.round() as i32);
            let prim_rotations_actual: Vec<_> = prim_rotations
                .iter()
                .map(|r| prim_trans_inv * r * point_group.prim_trans_mat)
                .collect();
            for rotation in prim_rotations_actual {
                assert!(prim_rotations_set.contains(&rotation));
            }

            // Compare rotations with the conventional
            let mut conv_rotations_set = HashSet::new();
            for rotation in traverse(&point_group_db.generators) {
                conv_rotations_set.insert(rotation);
            }
            let conv_trans_inv = point_group
                .conv_trans_mat
                .map(|e| e as f64)
                .try_inverse()
                .unwrap();
            let conv_rotations_actual: Vec<_> = prim_rotations
                .iter()
                .map(|r| {
                    let ret = conv_trans_inv
                        * r.map(|e| e as f64)
                        * point_group.conv_trans_mat.map(|e| e as f64);
                    ret.map(|e| e.round() as i32)
                })
                .collect();
            for rotation in conv_rotations_actual {
                assert!(conv_rotations_set.contains(&rotation));
            }
        }
    }
}
