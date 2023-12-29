use std::collections::HashMap;

use nalgebra::{Dyn, Matrix3, OMatrix, OVector, Vector3, U3};

use super::point_group::identify_point_group;
use crate::base::error::MoyoError;
use crate::base::operation::AbstractOperations;
use crate::base::tolerance::EPS;
use crate::base::transformation::{OriginShift, Transformation, TransformationMatrix};
use crate::data::arithmetic_crystal_class::{get_arithmetic_crystal_class_entry, ArithmeticNumber};
use crate::data::classification::CrystalSystem;
use crate::data::hall_symbol::HallSymbol;
use crate::data::hall_symbol_database::{get_hall_symbol_entry, HallNumber, Number};
use crate::data::point_group::PointGroupRepresentative;
use crate::identify::point_group::PointGroup;
use crate::math::snf::SNF;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Setting {
    HallNumber(HallNumber),
    /// The setting of the smallest Hall number
    Spglib,
    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral, and origin choice 2 for centrosymmetric space groups
    Standard,
}

impl Setting {
    pub fn hall_numbers(&self) -> Vec<HallNumber> {
        match self {
            Setting::HallNumber(hall_number) => vec![*hall_number],
            Setting::Spglib => vec![
                1, 2, 3, 6, 9, 18, 21, 30, 39, 57, 60, 63, 72, 81, 90, 108, 109, 112, 115, 116,
                119, 122, 123, 124, 125, 128, 134, 137, 143, 149, 155, 161, 164, 170, 173, 176,
                182, 185, 191, 197, 203, 209, 212, 215, 218, 221, 227, 228, 230, 233, 239, 245,
                251, 257, 263, 266, 269, 275, 278, 284, 290, 292, 298, 304, 310, 313, 316, 322,
                334, 335, 337, 338, 341, 343, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358,
                359, 361, 363, 364, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
                378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393,
                394, 395, 396, 397, 398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413,
                414, 416, 418, 419, 420, 422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436,
                438, 439, 440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456,
                457, 458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
                475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490,
                491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506, 507, 508,
                509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523, 524, 525, 527,
                529, 530,
            ],
            Setting::Standard => vec![
                1, 2, 3, 6, 9, 18, 21, 30, 39, 57, 60, 63, 72, 81, 90, 108, 109, 112, 115, 116,
                119, 122, 123, 124, 125, 128, 134, 137, 143, 149, 155, 161, 164, 170, 173, 176,
                182, 185, 191, 197, 203, 209, 212, 215, 218, 221, 227, 229, 230, 234, 239, 245,
                251, 257, 263, 266, 269, 275, 279, 284, 290, 292, 298, 304, 310, 313, 316, 323,
                334, 336, 337, 338, 341, 343, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358,
                360, 362, 363, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
                378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393,
                394, 395, 396, 397, 398, 399, 400, 401, 403, 405, 406, 407, 409, 411, 412, 413,
                415, 417, 418, 419, 421, 423, 424, 425, 427, 429, 430, 431, 432, 433, 435, 436,
                438, 439, 440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456,
                457, 458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
                475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490,
                491, 492, 493, 494, 496, 497, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508,
                509, 510, 511, 512, 513, 514, 515, 516, 517, 519, 520, 522, 523, 524, 526, 528,
                529, 530,
            ],
        }
    }
}

#[derive(Debug)]
pub struct SpaceGroup {
    pub number: Number,
    pub hall_number: HallNumber,
    /// Transformation to the representative for `hall_number` in primitive
    pub transformation: Transformation,
}

fn get_correction_transformation_matrices(
    arithmetic_number: ArithmeticNumber,
) -> Vec<TransformationMatrix> {
    let (_, _, geometric_crystal_class, _) = get_arithmetic_crystal_class_entry(arithmetic_number);
    let crystal_system = CrystalSystem::from_geometric_crystal_class(geometric_crystal_class);

    // conventional -> conventional(cell choice 1 for monoclinic, abc for orthorhombic)
    let convs = match crystal_system {
        CrystalSystem::Monoclinic => vec![
            TransformationMatrix::identity(),
            // b2 to b1
            TransformationMatrix::new(
                0, 0, -1, //
                0, 1, 0, //
                1, 0, -1, //
            ),
            // b3 to b1
            TransformationMatrix::new(
                -1, 0, 1, //
                0, 1, 0, //
                -1, 0, 0, //
            ),
        ],
        CrystalSystem::Orthorhombic => vec![
            // abc
            TransformationMatrix::identity(),
            // ba-c
            TransformationMatrix::new(
                0, 1, 0, //
                1, 0, 0, //
                0, 0, -1, //
            ),
            // cab
            TransformationMatrix::new(
                0, 0, 1, //
                1, 0, 0, //
                0, 1, 0, //
            ),
            // -cba
            TransformationMatrix::new(
                0, 0, -1, //
                0, 1, 0, //
                1, 0, 0, //
            ),
            // bca
            TransformationMatrix::new(
                0, 1, 0, //
                0, 0, 1, //
                1, 0, 0, //
            ),
            // a-cb
            TransformationMatrix::new(
                1, 0, 0, //
                0, 0, -1, //
                0, 1, 0, //
            ),
        ],
        _ => vec![TransformationMatrix::identity()],
    };

    // primitive -> conventional -> conventional(cell choice 1 for monoclinic, abc for orthorhombic) -> primitive
    let centering =
        PointGroupRepresentative::from_arithmetic_crystal_class(arithmetic_number).centering;
    let corrections: Vec<TransformationMatrix> = convs
        .iter()
        .map(|trans_corr| {
            let corr = centering.transformation_matrix().map(|e| e as f64)
                * trans_corr.map(|e| e as f64)
                * centering.inverse();
            corr.map(|e| e.round() as i32)
        })
        .filter(|corr| corr.map(|e| e as f64).determinant().round() as i32 == 1)
        .collect();
    corrections
}

pub fn identify_space_group(
    prim_operations: &AbstractOperations,
    setting: Setting,
) -> Result<SpaceGroup, MoyoError> {
    // point_group.trans_mat: self -> primitive
    let point_group = identify_point_group(&prim_operations.rotations)?;

    let new_prim_operations = prim_operations.transform(
        &point_group.prim_trans_mat.map(|e| e as f64),
        &OriginShift::zeros(),
    );
    let mut hm_translations = HashMap::new();
    for (rotation, translation) in new_prim_operations
        .rotations
        .iter()
        .zip(new_prim_operations.translations.iter())
    {
        hm_translations.insert(rotation.clone(), translation.clone());
    }
    dbg!(&hm_translations);

    for hall_number in setting.hall_numbers() {
        let entry = get_hall_symbol_entry(hall_number);
        if entry.arithmetic_number != point_group.arithmetic_number {
            continue;
        }

        let hall_symbol = HallSymbol::from_hall_number(hall_number);
        let other_prim_generators = hall_symbol.primitive_generators();

        // Solve (E, c)^-1 (R, t_target) (E, c) = (R, t_other) (mod 1) (for all (R, t_other) in other_prim_generators)
        // (R, R * c - c + t_target) = (R, t_other) (mod 1)
        // <-> (R - E) * c = t_other - t_target (mod 1)
        let mut a = OMatrix::<i32, Dyn, U3>::zeros(3 * other_prim_generators.rotations.len());
        let mut b = OVector::<f64, Dyn>::zeros(3 * other_prim_generators.rotations.len());
        for (k, (rotation, other_translation)) in other_prim_generators
            .rotations
            .iter()
            .zip(other_prim_generators.translations.iter())
            .enumerate()
        {
            let target_translation = hm_translations.get(rotation).unwrap();
            let ak = rotation - Matrix3::<i32>::identity();
            let bk = other_translation - target_translation;
            for i in 0..3 {
                for j in 0..3 {
                    a[(3 * k + i, j)] = ak[(i, j)];
                }
                b[3 * k + i] = bk[i];
            }
        }

        if let Some(origin_shift) = solve_mod1(&a, &b) {
            dbg!(hall_number);
            dbg!(a.map(|e| e as f64) * origin_shift.clone() - b);
            dbg!(&origin_shift);
            return Ok(SpaceGroup {
                number: entry.number,
                hall_number: hall_number,
                transformation: Transformation::new(point_group.prim_trans_mat, origin_shift),
            });
        }
    }

    Err(MoyoError::SpaceGroupTypeIdentificationError)
}

/// TODO: this function seems to be wrong!!!
/// Solve a * x = b (mod 1)
fn solve_mod1(a: &OMatrix<i32, Dyn, U3>, b: &OVector<f64, Dyn>) -> Option<Vector3<f64>> {
    // Solve snf.d * y = snf.l * b (x = snf.r * y)
    let snf = SNF::new(a);
    let mut y = Vector3::<f64>::zeros();
    let lb = snf.l.map(|e| e as f64) * b;
    for i in 0..3 {
        if snf.d[(i, i)] == 0 {
            let lbi = lb[i] - lb[i].round();
            if lbi.abs() > EPS {
                return None;
            }
        } else {
            y[i] = lb[i] / (snf.d[(i, i)] as f64);
        }
    }
    let mut x = snf.r.map(|e| e as f64) * y;
    x -= x.map(|e| e.round()); // mod 1
    Some(x)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use nalgebra::{matrix, vector, Dyn, OMatrix, OVector, RowVector3, U3};

    use crate::base::transformation::OriginShift;
    use crate::data::hall_symbol::HallSymbol;
    use crate::data::hall_symbol_database::get_hall_symbol_entry;

    use super::{
        get_correction_transformation_matrices, identify_space_group, solve_mod1, Setting,
    };

    #[test]
    fn test_solve_mod1() {
        let a = OMatrix::<i32, Dyn, U3>::from_rows(&[
            RowVector3::new(1, 0, 0),
            RowVector3::new(0, 1, 0),
            RowVector3::new(0, 0, 1),
        ]);
        let b = OVector::<f64, Dyn>::from_row_slice(&[1.0, 0.0, 0.0]);
        let x = solve_mod1(&a, &b);
        // TODO
    }

    #[test]
    fn test_get_correction_transformation_matrices() {
        let hall_number = 21; // P 1 c 1
        let hall_symbol = HallSymbol::from_hall_number(hall_number);
        let operations = hall_symbol.traverse();

        // conventional -> primitive
        let linear = hall_symbol.lattice_symbol.inverse();
        let prim_operations = operations.transform(&linear, &OriginShift::zeros());

        // The correction transformation matrices should change the group into P1c1, P1a1, and P1n1
        let entry = get_hall_symbol_entry(hall_number);
        let corrections = get_correction_transformation_matrices(entry.arithmetic_number);
        let expects = vec![
            vector![0.0, 0.0, 0.5],
            vector![0.5, 0.0, 0.0],
            vector![-0.5, 0.0, -0.5],
        ];
        for (i, corr) in corrections.iter().enumerate() {
            let corr_prim_operations =
                prim_operations.transform(&corr.map(|e| e as f64), &OriginShift::zeros());
            let mut hm_translations = HashMap::new();
            for (rotation, translation) in corr_prim_operations
                .rotations
                .iter()
                .zip(corr_prim_operations.translations.iter())
            {
                hm_translations.insert(rotation.clone(), translation.clone());
            }
            let r = matrix![
                1, 0, 0;
                0, -1, 0;
                0, 0, 1;
            ];
            assert_relative_eq!(*hm_translations.get(&r).unwrap(), expects[i]);
        }
    }

    #[test]
    fn test_identify_space_group() {
        // for hall_number in 1..=530 {
        for hall_number in [21] {
            dbg!(hall_number);
            let hall_symbol = HallSymbol::from_hall_number(hall_number);
            let operations = hall_symbol.traverse();

            // conventional -> primitive
            let linear = hall_symbol.lattice_symbol.inverse();
            let prim_operations = operations.transform(&linear, &OriginShift::zeros());

            let space_group = identify_space_group(&prim_operations, Setting::Spglib).unwrap();

            let entry = get_hall_symbol_entry(hall_number);
            assert_eq!(space_group.number, entry.number);
        }
    }
}
