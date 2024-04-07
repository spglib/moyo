use std::cmp::Ordering;

use itertools::Itertools;
use log::debug;
use nalgebra::{Dyn, Matrix3, OMatrix, OVector, U9};

use crate::base::{MoyoError, Rotation, UnimodularLinear};
use crate::data::{
    iter_arithmetic_crystal_entry, ArithmeticNumber, Centering, CrystalSystem,
    GeometricCrystalClass, PointGroupRepresentative,
};
use crate::math::IntegerLinearSystem;

/// Crystallographic point group with group-type information
#[derive(Debug)]
pub struct PointGroup {
    pub arithmetic_number: ArithmeticNumber,
    /// Transformation matrix to the representative for `arithmetic_number` in the primitive basis
    pub prim_trans_mat: UnimodularLinear,
}

impl PointGroup {
    /// Identify the arithmetic crystal class from the given rotations and transformation matrix into the representative
    /// Assume the rotations are given in the (reduced) primitive basis
    pub fn new(prim_rotations: &Vec<Rotation>) -> Result<Self, MoyoError> {
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

            for (arithmetic_crystal_class, point_group_db) in
                arithmetic_crystal_class_candidates.iter()
            {
                let centering = point_group_db.centering;
                if centering.order() as i32 != det {
                    continue;
                }
                // conv_trans_mat: self -> conventional
                let prim_trans_mat =
                    (conv_trans_mat * centering.inverse()).map(|e| e.round() as i32);
                if prim_trans_mat.map(|e| e as f64).determinant().round() as i32 != 1 {
                    return Err(MoyoError::ArithmeticCrystalClassIdentificationError);
                }

                return Ok(PointGroup {
                    arithmetic_number: *arithmetic_crystal_class,
                    prim_trans_mat,
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

    for entry in iter_arithmetic_crystal_entry() {
        if entry.geometric_crystal_class != geometric_crystal_class {
            continue;
        }

        let point_group_db =
            PointGroupRepresentative::from_arithmetic_crystal_class(entry.arithmetic_number);
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
                    let mut prim_trans_mat = UnimodularLinear::zeros();
                    for (i, matrix) in trans_mat_basis.iter().enumerate() {
                        prim_trans_mat += comb[i] * matrix;
                    }
                    let det = prim_trans_mat.map(|e| e as f64).determinant().round() as i32;
                    if det == 1 {
                        return Ok(PointGroup {
                            arithmetic_number: entry.arithmetic_number,
                            prim_trans_mat,
                        });
                    }
                }
            }
        }
    }

    Err(MoyoError::ArithmeticCrystalClassIdentificationError)
}

#[allow(dead_code)]
pub fn integral_normalizer(
    prim_rotations: &Vec<Rotation>,
    prim_generators: &Vec<Rotation>,
) -> Vec<UnimodularLinear> {
    let rotation_types = prim_rotations
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();
    let prim_generator_rotation_types = prim_generators
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();

    // Try to map generators
    let order = prim_rotations.len();
    let candidates: Vec<Vec<usize>> = prim_generator_rotation_types
        .iter()
        .map(|&rotation_type| {
            (0..order)
                .filter(|&i| rotation_types[i] == rotation_type)
                .collect()
        })
        .collect();

    // TODO: unify with match_with_point_group
    let mut conjugators = vec![];
    for pivot in candidates
        .iter()
        .map(|v| v.iter())
        .multi_cartesian_product()
    {
        // Solve P^-1 * prim_rotations[prim_rotations[pivot[i]]] * P = prim_generators[i] (for all i)
        if let Some(trans_mat_basis) = sylvester(
            &pivot
                .iter()
                .map(|&i| prim_rotations[*i])
                .collect::<Vec<_>>(),
            &prim_generators,
        ) {
            // Search integer linear combination such that the transformation matrix is unimodular
            // Consider coefficients in [-2, 2], which will be sufficient for Delaunay reduced basis
            for comb in (0..trans_mat_basis.len())
                .map(|_| -2..=2)
                .multi_cartesian_product()
            {
                let mut prim_trans_mat = UnimodularLinear::zeros();
                for (i, matrix) in trans_mat_basis.iter().enumerate() {
                    prim_trans_mat += comb[i] * matrix;
                }
                let det = prim_trans_mat.map(|e| e as f64).determinant().round() as i32;
                if det < 0 {
                    prim_trans_mat *= -1;
                }
                if det == 1 {
                    conjugators.push(prim_trans_mat);
                    break;
                }
            }
        }
    }
    conjugators
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
                UnimodularLinear::new(
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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::{integral_normalizer, PointGroup, PointGroupRepresentative};
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

    #[test]
    fn test_integral_normalizer() {
        for arithmetic_number in 1..=73 {
            let point_group_db =
                PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
            let prim_generators = point_group_db.primitive_generators();
            let prim_rotations = traverse(&prim_generators);

            let mut prim_rotations_set = HashSet::new();
            for rotation in prim_rotations.iter() {
                prim_rotations_set.insert(rotation.clone());
            }

            let conjugators = integral_normalizer(&prim_rotations, &prim_generators);
            assert!(conjugators.len() > 0);

            for linear in conjugators.iter() {
                let linear_inv = linear
                    .map(|e| e as f64)
                    .try_inverse()
                    .unwrap()
                    .map(|e| e.round() as i32);
                let prim_rotations_actual: Vec<_> = prim_rotations
                    .iter()
                    .map(|r| linear_inv * r * linear)
                    .collect();
                for rotation in prim_rotations_actual {
                    assert!(prim_rotations_set.contains(&rotation));
                }
            }
        }
    }
}
