use std::cmp::Ordering;

use itertools::Itertools;
use log::debug;
use nalgebra::Matrix3;
use serde::{Deserialize, Serialize};

use super::rotation_type::{identify_rotation_type, RotationType};
use crate::base::{
    project_rotations, Lattice, MoyoError, Operation, Rotations, Translation, UnimodularLinear,
    UnimodularTransformation,
};
use crate::data::{
    iter_arithmetic_crystal_entry, ArithmeticNumber, Centering, CrystalSystem,
    GeometricCrystalClass, PointGroupRepresentative,
};
use crate::math::sylvester3;

/// Crystallographic point group with group-type information
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PointGroup {
    pub arithmetic_number: ArithmeticNumber,
    /// Transformation matrix to the representative for `arithmetic_number` in the primitive basis
    pub prim_trans_mat: UnimodularLinear,
}

impl PointGroup {
    /// Given rotations, identify the arithmetic crystal class and transformation matrix into the representative
    /// Assume the rotations are given in the (reduced) primitive basis
    pub fn new(prim_rotations: &Rotations) -> Result<Self, MoyoError> {
        let rotation_types = prim_rotations.iter().map(identify_rotation_type).collect();
        let geometric_crystal_class = identify_geometric_crystal_class(&rotation_types)?;
        debug!("Geometric crystal class: {:?}", geometric_crystal_class);

        let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);
        match crystal_system {
            // Skip triclinic cases as trivial
            CrystalSystem::Triclinic => match geometric_crystal_class {
                GeometricCrystalClass::C1 => Ok(PointGroup {
                    arithmetic_number: 1,
                    prim_trans_mat: UnimodularLinear::identity(),
                }),
                GeometricCrystalClass::Ci => Ok(PointGroup {
                    arithmetic_number: 2,
                    prim_trans_mat: UnimodularLinear::identity(),
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

    /// Given lattice and rotations, identify the arithmetic crystal class and transformation matrix into the representative
    pub fn from_lattice(lattice: &Lattice, prim_rotations: &Rotations) -> Result<Self, MoyoError> {
        let (_, reduced_trans_mat) = lattice.minkowski_reduce()?;
        let reduced_prim_operations = UnimodularTransformation::from_linear(reduced_trans_mat)
            .transform_operations(
                &prim_rotations
                    .iter()
                    .map(|r| Operation::new(*r, Translation::zeros()))
                    .collect::<Vec<_>>(),
            );
        let reduced_prim_rotations = project_rotations(&reduced_prim_operations);

        let reduced_point_group = Self::new(&reduced_prim_rotations)?;
        Ok(PointGroup {
            arithmetic_number: reduced_point_group.arithmetic_number,
            prim_trans_mat: reduced_point_group.prim_trans_mat * reduced_trans_mat,
        })
    }
}

/// Faster matching algorithm for cubic point groups
fn match_with_cubic_point_group(
    prim_rotations: &Rotations,
    rotation_types: &[RotationType],
    geometric_crystal_class: GeometricCrystalClass,
) -> Result<PointGroup, MoyoError> {
    let arithmetic_crystal_class_candidates = iter_arithmetic_crystal_entry()
        .filter_map(|entry| {
            if entry.geometric_crystal_class != geometric_crystal_class {
                return None;
            }
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(entry.arithmetic_number);
            Some((entry.arithmetic_number, point_group_db))
        })
        .collect::<Vec<_>>();
    let primitive_arithmetic_crystal_class = arithmetic_crystal_class_candidates
        .iter()
        .find(|(_, point_group_db)| point_group_db.centering == Centering::P)
        .unwrap();
    let other_prim_generators = primitive_arithmetic_crystal_class.1.primitive_generators();

    for trans_mat_basis in iter_trans_mat_basis(
        prim_rotations.clone(),
        rotation_types.to_vec(),
        other_prim_generators,
    ) {
        // conv_trans_mat: self -> P-centering (primitive)
        // The dimension of linear integer system should be one for cubic.
        assert_eq!(trans_mat_basis.len(), 1);
        let mut conv_trans_mat = trans_mat_basis[0].map(|e| e as f64);

        // Guarantee det > 0
        let mut det = conv_trans_mat.determinant().round() as i32;
        match det.cmp(&0) {
            Ordering::Less => {
                conv_trans_mat *= -1.0;
                det *= -1;
            }
            Ordering::Equal => continue,
            Ordering::Greater => {}
        }

        for (arithmetic_crystal_class, point_group_db) in arithmetic_crystal_class_candidates.iter()
        {
            let centering = point_group_db.centering;
            if centering.order() as i32 != det {
                continue;
            }
            // conv_trans_mat: self -> conventional
            let prim_trans_mat = (conv_trans_mat * centering.inverse()).map(|e| e.round() as i32);
            if prim_trans_mat.map(|e| e as f64).determinant().round() as i32 != 1 {
                return Err(MoyoError::ArithmeticCrystalClassIdentificationError);
            }

            return Ok(PointGroup {
                arithmetic_number: *arithmetic_crystal_class,
                prim_trans_mat,
            });
        }
    }

    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

fn match_with_point_group(
    prim_rotations: &Rotations,
    rotation_types: &[RotationType],
    geometric_crystal_class: GeometricCrystalClass,
) -> Result<PointGroup, MoyoError> {
    for entry in iter_arithmetic_crystal_entry() {
        if entry.geometric_crystal_class != geometric_crystal_class {
            continue;
        }

        let point_group_db =
            PointGroupRepresentative::from_arithmetic_crystal_class(entry.arithmetic_number);
        let other_prim_generators = point_group_db.primitive_generators();

        for trans_mat_basis in iter_trans_mat_basis(
            prim_rotations.clone(),
            rotation_types.to_vec(),
            other_prim_generators,
        ) {
            if let Some(prim_trans_mat) = iter_unimodular_trans_mat(trans_mat_basis).nth(0) {
                return Ok(PointGroup {
                    arithmetic_number: entry.arithmetic_number,
                    prim_trans_mat,
                });
            }
        }
    }

    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

/// Use look up table in Table 6 of https://arxiv.org/pdf/1808.01590.pdf
fn identify_geometric_crystal_class(
    rotation_types: &Vec<RotationType>,
) -> Result<GeometricCrystalClass, MoyoError> {
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
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0] => Ok(GeometricCrystalClass::C1),
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0] => Ok(GeometricCrystalClass::Ci),
        // Monoclinic
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0] => Ok(GeometricCrystalClass::C2),
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0] => Ok(GeometricCrystalClass::C1h),
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0] => Ok(GeometricCrystalClass::C2h),
        // Orthorhombic
        [0, 0, 0, 0, 0, 1, 3, 0, 0, 0] => Ok(GeometricCrystalClass::D2),
        [0, 0, 0, 2, 0, 1, 1, 0, 0, 0] => Ok(GeometricCrystalClass::C2v),
        [0, 0, 0, 3, 1, 1, 3, 0, 0, 0] => Ok(GeometricCrystalClass::D2h),
        // Tetragonal
        [0, 0, 0, 0, 0, 1, 1, 0, 2, 0] => Ok(GeometricCrystalClass::C4),
        [0, 2, 0, 0, 0, 1, 1, 0, 0, 0] => Ok(GeometricCrystalClass::S4),
        [0, 2, 0, 1, 1, 1, 1, 0, 2, 0] => Ok(GeometricCrystalClass::C4h),
        [0, 0, 0, 0, 0, 1, 5, 0, 2, 0] => Ok(GeometricCrystalClass::D4),
        [0, 0, 0, 4, 0, 1, 1, 0, 2, 0] => Ok(GeometricCrystalClass::C4v),
        [0, 2, 0, 2, 0, 1, 3, 0, 0, 0] => Ok(GeometricCrystalClass::D2d),
        [0, 2, 0, 5, 1, 1, 5, 0, 2, 0] => Ok(GeometricCrystalClass::D4h),
        // Trigonal
        [0, 0, 0, 0, 0, 1, 0, 2, 0, 0] => Ok(GeometricCrystalClass::C3),
        [0, 0, 2, 0, 1, 1, 0, 2, 0, 0] => Ok(GeometricCrystalClass::C3i),
        [0, 0, 0, 0, 0, 1, 3, 2, 0, 0] => Ok(GeometricCrystalClass::D3),
        [0, 0, 0, 3, 0, 1, 0, 2, 0, 0] => Ok(GeometricCrystalClass::C3v),
        [0, 0, 2, 3, 1, 1, 3, 2, 0, 0] => Ok(GeometricCrystalClass::D3d),
        // Hexagonal
        [0, 0, 0, 0, 0, 1, 1, 2, 0, 2] => Ok(GeometricCrystalClass::C6),
        [2, 0, 0, 1, 0, 1, 0, 2, 0, 0] => Ok(GeometricCrystalClass::C3h),
        [2, 0, 2, 1, 1, 1, 1, 2, 0, 2] => Ok(GeometricCrystalClass::C6h),
        [0, 0, 0, 0, 0, 1, 7, 2, 0, 2] => Ok(GeometricCrystalClass::D6),
        [0, 0, 0, 6, 0, 1, 1, 2, 0, 2] => Ok(GeometricCrystalClass::C6v),
        [2, 0, 0, 4, 0, 1, 3, 2, 0, 0] => Ok(GeometricCrystalClass::D3h),
        [2, 0, 2, 7, 1, 1, 7, 2, 0, 2] => Ok(GeometricCrystalClass::D6h),
        // Cubic
        [0, 0, 0, 0, 0, 1, 3, 8, 0, 0] => Ok(GeometricCrystalClass::T),
        [0, 0, 8, 3, 1, 1, 3, 8, 0, 0] => Ok(GeometricCrystalClass::Th),
        [0, 0, 0, 0, 0, 1, 9, 8, 6, 0] => Ok(GeometricCrystalClass::O),
        [0, 6, 0, 6, 0, 1, 3, 8, 0, 0] => Ok(GeometricCrystalClass::Td),
        [0, 6, 8, 9, 1, 1, 9, 8, 6, 0] => Ok(GeometricCrystalClass::Oh),
        // Unknown
        _ => {
            debug!(
                "Unknown geometric crystal class: {:?}",
                rotation_types_count
            );
            Err(MoyoError::GeometricCrystalClassIdentificationError)
        }
    }
}

/// Return iterator of basis for transformation matrices that map the point group `prim_rotations` to a point group formed by `other_prim_rotation_generators`.
pub fn iter_trans_mat_basis(
    prim_rotations: Rotations,
    rotation_types: Vec<RotationType>,
    other_prim_rotation_generators: Rotations,
) -> impl Iterator<Item = Vec<Matrix3<i32>>> {
    let other_prim_generator_rotation_types = other_prim_rotation_generators
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();

    let order = prim_rotations.len();
    let candidates: Vec<Vec<usize>> = other_prim_generator_rotation_types
        .iter()
        .map(|&rotation_type| {
            (0..order)
                .filter(|&i| rotation_types[i] == rotation_type)
                .collect()
        })
        .collect();

    candidates
        .into_iter()
        .map(|v| v.into_iter())
        .multi_cartesian_product()
        .filter_map(move |pivot| {
            // Solve P^-1 * prim_rotations[pivot[i]] * P = prim_rotation_generators[i] (for all i)
            sylvester3(
                &pivot.iter().map(|&i| prim_rotations[i]).collect::<Vec<_>>(),
                &other_prim_rotation_generators,
            )
        })
}

/// Search integer linear combination such that the transformation matrix is unimodular
/// Consider coefficients in [-2, 2], which will be sufficient for Delaunay reduced basis
pub fn iter_unimodular_trans_mat(
    trans_mat_basis: Vec<Matrix3<i32>>,
) -> impl Iterator<Item = UnimodularLinear> {
    // First try with coefficients in [-1, 1]
    let iter_multi_1 = (0..trans_mat_basis.len())
        .map(|_| -1..=1)
        .multi_cartesian_product();
    let iter_multi_2 = (0..trans_mat_basis.len())
        .map(|_| -2..=2)
        .multi_cartesian_product()
        .filter(|comb| comb.iter().any(|&e| (e as i32).abs() == 2));

    iter_multi_1.chain(iter_multi_2).filter_map(move |comb| {
        // prim_trans_mat: self -> DB(primitive)
        let mut prim_trans_mat = UnimodularLinear::zeros();
        for (i, matrix) in trans_mat_basis.iter().enumerate() {
            prim_trans_mat += comb[i] * matrix;
        }
        let det = prim_trans_mat.map(|e| e as f64).determinant().round() as i32;
        if det == 1 {
            Some(prim_trans_mat)
        } else {
            None
        }
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::base::traverse;

    #[test]
    fn test_point_group_match() {
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
        }
    }
}
