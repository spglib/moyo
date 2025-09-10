use std::collections::HashMap;

use log::debug;
use nalgebra::{Dyn, Matrix3, OMatrix, OVector, Vector3, U3};

use super::point_group::PointGroup;
use crate::base::{
    project_rotations, Lattice, MoyoError, Operations, OriginShift, UnimodularLinear,
    UnimodularTransformation,
};
use crate::data::{
    arithmetic_crystal_class_entry, hall_symbol_entry, ArithmeticNumber, GeometricCrystalClass,
    HallNumber, HallSymbol, Number, PointGroupRepresentative, Setting,
};
use crate::math::SNF;

#[derive(Debug, Clone)]
pub struct SpaceGroup {
    pub number: Number,
    pub hall_number: HallNumber,
    /// Transformation to the representative for `hall_number` in primitive
    pub transformation: UnimodularTransformation,
}

impl SpaceGroup {
    /// Identify the space group from the primitive operations.
    /// epsilon: tolerance for comparing translation parts
    pub fn new(
        prim_operations: &Operations,
        setting: Setting,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        // point_group.trans_mat: self -> primitive
        let prim_rotations = project_rotations(prim_operations);
        let point_group = PointGroup::new(&prim_rotations)?;
        debug!(
            "Arithmetic crystal class: No. {}",
            point_group.arithmetic_number
        );

        for hall_number in setting.hall_numbers() {
            let entry = hall_symbol_entry(hall_number).unwrap();
            if entry.arithmetic_number != point_group.arithmetic_number {
                continue;
            }

            let hall_symbol = HallSymbol::from_hall_number(hall_number)
                .ok_or(MoyoError::SpaceGroupTypeIdentificationError)?;
            let db_prim_generators = hall_symbol.primitive_generators();

            // Try several correction transformation matrices for monoclinic and orthorhombic
            for trans_mat_corr in correction_transformation_matrices(entry.arithmetic_number) {
                let trans_mat = point_group.prim_trans_mat * trans_mat_corr;
                if let Some(origin_shift) =
                    match_origin_shift(prim_operations, &trans_mat, &db_prim_generators, epsilon)
                {
                    debug!(
                        "Matched with Hall number {} (No. {})",
                        hall_number, entry.number
                    );
                    return Ok(Self {
                        number: entry.number,
                        hall_number,
                        transformation: UnimodularTransformation::new(trans_mat, origin_shift),
                    });
                }
            }
        }

        Err(MoyoError::SpaceGroupTypeIdentificationError)
    }

    pub fn from_lattice(
        lattice: &Lattice,
        prim_operations: &Operations,
        setting: Setting,
        epsilon: f64,
    ) -> Result<Self, MoyoError> {
        let (_, reduced_trans_mat) = lattice.minkowski_reduce()?;
        let to_reduced = UnimodularTransformation::from_linear(reduced_trans_mat);
        let reduced_prim_operations = to_reduced.transform_operations(prim_operations);

        let reduced_space_group = Self::new(&reduced_prim_operations, setting, epsilon)?;
        Ok(SpaceGroup {
            number: reduced_space_group.number,
            hall_number: reduced_space_group.hall_number,
            transformation: reduced_space_group.transformation * to_reduced,
        })
    }

    pub fn from_hall_number_and_transformation(
        hall_number: HallNumber,
        transformation: UnimodularTransformation,
    ) -> Result<Self, MoyoError> {
        let entry =
            hall_symbol_entry(hall_number).ok_or(MoyoError::SpaceGroupTypeIdentificationError)?;
        Ok(Self {
            number: entry.number,
            hall_number,
            transformation,
        })
    }
}

fn correction_transformation_matrices(
    arithmetic_number: ArithmeticNumber,
) -> Vec<UnimodularLinear> {
    let geometric_crystal_class = arithmetic_crystal_class_entry(arithmetic_number)
        .unwrap()
        .geometric_crystal_class;

    // conventional -> conventional(standard)
    let convs = match geometric_crystal_class {
        // Monoclinic crystal system
        GeometricCrystalClass::C2 | GeometricCrystalClass::C1h | GeometricCrystalClass::C2h => {
            vec![
                UnimodularLinear::identity(),
                // b2 to b1
                UnimodularLinear::new(
                    0, 0, -1, //
                    0, 1, 0, //
                    1, 0, -1, //
                ),
                // b3 to b1
                UnimodularLinear::new(
                    -1, 0, 1, //
                    0, 1, 0, //
                    -1, 0, 0, //
                ),
            ]
        }
        // Orthorhombic crystal system
        GeometricCrystalClass::D2 | GeometricCrystalClass::C2v | GeometricCrystalClass::D2h => {
            vec![
                // abc
                UnimodularLinear::identity(),
                // ba-c
                UnimodularLinear::new(
                    0, 1, 0, //
                    1, 0, 0, //
                    0, 0, -1, //
                ),
                // cab
                UnimodularLinear::new(
                    0, 0, 1, //
                    1, 0, 0, //
                    0, 1, 0, //
                ),
                // -cba
                UnimodularLinear::new(
                    0, 0, -1, //
                    0, 1, 0, //
                    1, 0, 0, //
                ),
                // bca
                UnimodularLinear::new(
                    0, 1, 0, //
                    0, 0, 1, //
                    1, 0, 0, //
                ),
                // a-cb
                UnimodularLinear::new(
                    1, 0, 0, //
                    0, 0, -1, //
                    0, 1, 0, //
                ),
            ]
        }
        // m-3
        GeometricCrystalClass::Th => {
            vec![
                UnimodularLinear::identity(),
                UnimodularLinear::new(
                    0, 0, 1, //
                    0, -1, 0, //
                    1, 0, 0, //
                ),
            ]
        }
        _ => vec![UnimodularLinear::identity()],
    };

    // primitive -> conventional -> conventional(standard) -> primitive
    let point_group = PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number);
    let centering = point_group.centering;
    let corrections: Vec<UnimodularLinear> = convs
        .iter()
        .map(|trans_corr| {
            let corr = (centering.linear() * trans_corr).map(|e| e as f64) * centering.inverse();
            corr.map(|e| e.round() as i32)
        })
        .filter(|corr| corr.map(|e| e as f64).determinant().round() as i32 == 1)
        .collect();
    corrections
}

/// Search for origin_shift such that (trans_mat, origin_shift) transforms <prim_operations> into <db_prim_generators>
pub fn match_origin_shift(
    prim_operations: &Operations,
    trans_mat: &UnimodularLinear,
    db_prim_generators: &Operations,
    epsilon: f64,
) -> Option<OriginShift> {
    let new_prim_operations =
        UnimodularTransformation::from_linear(*trans_mat).transform_operations(prim_operations);
    let mut hm_translations = HashMap::new();
    for operation in new_prim_operations.iter() {
        hm_translations.insert(operation.rotation, operation.translation);
    }

    // Find origin_shift `c`: (P, c)^-1 G (P, c) = G_db
    //     (P, c) = (P, 0) (P, 0)^-1 (P, c) = (P, 0) (E, P^-1 c)
    //     s := P^-1 c
    //     G' := (P, 0)^-1 G (P, 0)
    //     (E, s)^-1 G' (E, s) = G_db
    // Solve (E, s)^-1 (R, t_target) (E, s) = (R, t_db) (mod 1) (for all (R, t_db) in db_prim_generators)
    //     (R, R * s - s + t_target) = (R, t_db) (mod 1)
    //     <-> (R - E) * s = t_db - t_target (mod 1)
    let mut a = OMatrix::<i32, Dyn, U3>::zeros(3 * db_prim_generators.len());
    let mut b = OVector::<f64, Dyn>::zeros(3 * db_prim_generators.len());
    for (k, operation) in db_prim_generators.iter().enumerate() {
        // Correction transformation matrix may not be normalizer of the point group. For example, mm2 -> 2mm
        let target_translation = hm_translations.get(&operation.rotation)?;

        let ak = operation.rotation - Matrix3::<i32>::identity();
        let bk = operation.translation - target_translation;
        for i in 0..3 {
            for j in 0..3 {
                a[(3 * k + i, j)] = ak[(i, j)];
            }
            b[3 * k + i] = bk[i];
        }
    }

    match solve_mod1(&a, &b, epsilon) {
        Some(s) => {
            let origin_shift = (trans_mat.map(|e| e as f64) * s).map(|e| e % 1.);
            Some(origin_shift)
        }
        None => None,
    }
}

/// Solve a * x = b (mod 1)
pub fn solve_mod1(
    a: &OMatrix<i32, Dyn, U3>,
    b: &OVector<f64, Dyn>,
    epsilon: f64,
) -> Option<Vector3<f64>> {
    // Solve snf.d * y = snf.l * b (x = snf.r * y)
    let snf = SNF::new(a);
    let mut y = Vector3::<f64>::zeros();
    let lb = snf.l.map(|e| e as f64) * b;
    for i in 0..3 {
        if snf.d[(i, i)] == 0 {
            let lbi = lb[i] - lb[i].round();
            if lbi.abs() > epsilon {
                return None;
            }
        } else {
            y[i] = lb[i] / (snf.d[(i, i)] as f64);
        }
    }

    let x = (snf.r.map(|e| e as f64) * y).map(|e| e % 1.);

    // Check solution
    let mut residual = a.map(|e| e as f64) * x - b;
    residual -= residual.map(|e| e.round()); // in [-0.5, 0.5]
    if residual.iter().any(|e| e.abs() > epsilon) {
        return None;
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use nalgebra::{matrix, vector, Dyn, OMatrix, OVector, RowVector3, U3};
    use rstest::rstest;
    use std::collections::HashMap;

    use crate::base::{UnimodularTransformation, EPS};
    use crate::data::{hall_symbol_entry, HallSymbol, Setting};

    use super::{correction_transformation_matrices, solve_mod1, SpaceGroup};

    #[test]
    fn test_solve_mod1() {
        let a = OMatrix::<i32, Dyn, U3>::from_rows(&[
            RowVector3::new(-2, 0, 0),
            RowVector3::new(0, -2, 0),
            RowVector3::new(0, 0, -2),
            RowVector3::new(-2, 0, 0),
            RowVector3::new(0, 0, 0),
            RowVector3::new(0, 0, -2),
        ]);
        let b = OVector::<f64, Dyn>::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.5, 0.0]);
        let x = solve_mod1(&a, &b, EPS);
        assert_eq!(x, None);
    }

    #[test]
    fn test_correction_transformation_matrices() {
        let hall_number = 21; // P 1 c 1
        let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        let prim_operations = hall_symbol.primitive_traverse();

        // The correction transformation matrices should change the group into P1c1, P1a1, and P1n1
        let entry = hall_symbol_entry(hall_number).unwrap();
        let corrections = correction_transformation_matrices(entry.arithmetic_number);
        let expects = vec![
            vector![0.0, 0.0, 0.5],
            vector![0.5, 0.0, 0.0],
            vector![-0.5, 0.0, -0.5],
        ];
        for (i, corr) in corrections.iter().enumerate() {
            let corr_prim_operations =
                UnimodularTransformation::from_linear(*corr).transform_operations(&prim_operations);

            let mut hm_translations = HashMap::new();
            for operation in corr_prim_operations.iter() {
                hm_translations.insert(operation.rotation, operation.translation);
            }
            let r = matrix![
                1, 0, 0;
                0, -1, 0;
                0, 0, 1;
            ];
            assert_relative_eq!(*hm_translations.get(&r).unwrap(), expects[i]);
        }
    }

    #[rstest]
    #[case(Setting::Spglib)]
    #[case(Setting::Standard)]
    fn test_identify_space_group(#[case] setting: Setting) {
        for hall_number in 1..=530 {
            let hall_symbol = HallSymbol::from_hall_number(hall_number).unwrap();
            let prim_operations = hall_symbol.primitive_traverse();

            let space_group = SpaceGroup::new(&prim_operations, setting, 1e-8).unwrap();

            // Check space group type
            let entry = hall_symbol_entry(hall_number).unwrap();
            assert_eq!(space_group.number, entry.number);

            // Check transformation matrix
            assert_eq!(
                space_group
                    .transformation
                    .linear_as_f64()
                    .determinant()
                    .round() as i32,
                1
            );

            let matched_hall_symbol =
                HallSymbol::from_hall_number(space_group.hall_number).unwrap();
            let matched_prim_operations = matched_hall_symbol.primitive_traverse();

            let mut hm_translations = HashMap::new();
            for operation in matched_prim_operations.iter() {
                hm_translations.insert(operation.rotation, operation.translation);
            }

            // Check transformation
            let transformed_prim_operations = space_group
                .transformation
                .transform_operations(&prim_operations);
            assert_eq!(
                matched_prim_operations.len(),
                transformed_prim_operations.len()
            );
            for operation in transformed_prim_operations.iter() {
                assert!(hm_translations.contains_key(&operation.rotation));
                let mut diff =
                    *hm_translations.get(&operation.rotation).unwrap() - operation.translation;
                diff -= diff.map(|e| e.round()); // in [-0.5, 0.5]
                assert_relative_eq!(diff, vector![0.0, 0.0, 0.0])
            }
        }
    }
}
