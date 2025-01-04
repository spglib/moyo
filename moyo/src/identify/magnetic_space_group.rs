use std::collections::HashMap;

use log::debug;
use nalgebra::{Dyn, OMatrix, OVector, U3};

use super::space_group::{solve_mod1, SpaceGroup};
use crate::base::{
    MagneticOperations, MoyoError, Operations, Rotation, Translation, UnimodularLinear,
    UnimodularTransformation,
};
use crate::data::{
    get_magnetic_space_group_type, uni_number_range, ConstructType, MagneticHallSymbol, Setting,
    UNINumber,
};

#[derive(Debug)]
pub struct MagneticSpaceGroup {
    pub uni_number: UNINumber,
    /// Transformation to the representative for `uni_number` in primitive
    pub transformation: UnimodularTransformation,
}

impl MagneticSpaceGroup {
    /// Identify the magnetic space group type from the primitive magnetic operations.
    /// epsilon: tolerance for comparing translation parts
    pub fn new(prim_mag_operations: &MagneticOperations, epsilon: f64) -> Result<Self, MoyoError> {
        let (ref_spg, construct_type) =
            identify_reference_space_group(prim_mag_operations, epsilon)
                .ok_or(MoyoError::ConstructTypeIdentificationError)?;
        debug!("Construct type: {:?}", construct_type);
        let setting = Setting::Standard;
        // std_ref_spg.transformation: primitive input -> primitive BNS setting
        let std_ref_spg = SpaceGroup::new(&ref_spg, setting, epsilon)?;
        debug!("Reference space group: {:?}", std_ref_spg.number);

        let uni_number_range = uni_number_range(std_ref_spg.number)
            .ok_or(MoyoError::SpaceGroupTypeIdentificationError)?;
        for uni_number in uni_number_range {
            if get_magnetic_space_group_type(uni_number)
                .unwrap()
                .construct_type
                != construct_type
            {
                continue;
            }

            if construct_type == ConstructType::Type1 || construct_type == ConstructType::Type2 {
                // No need to further check the magnetic operations
                return Ok(Self {
                    uni_number,
                    transformation: std_ref_spg.transformation,
                });
            }

            // TODO:

            let mhs = MagneticHallSymbol::from_uni_number(uni_number)
                .ok_or(MoyoError::MagneticSpaceGroupTypeIdentificationError)?;
            let db_prim_mag_generators = mhs.primitive_generators();
            let db_ref_prim_generators =
                Self::get_db_reference_space_group_primitive_generators(&mhs, &construct_type);

            // The correction transformations keep the reference space group of `tmp_prim_mag_operations`
            // TODO: precompute the correction transformation matrices
            let correction_transformations = todo!();
            // for corr_trans_mat in correction_transformation_matrices {
            //     // (trans_mat, origin_shift): primitive input -> primitive DB
            //     let trans_mat = std_ref_spg.transformation.linear * corr_trans_mat;
            //     if let Some(origin_shift) = Self::match_origin_shift(
            //         prim_mag_operations,
            //         &trans_mat,
            //         &db_prim_mag_generators,
            //         epsilon,
            //     ) {
            //         debug!("Matched with UNI number {}", uni_number);
            //         return Ok(Self {
            //             uni_number,
            //             transformation: UnimodularTransformation::new(trans_mat, origin_shift),
            //         });
            //     }
            // }
        }
        Err(MoyoError::MagneticSpaceGroupTypeIdentificationError)
    }

    fn get_db_reference_space_group_primitive_generators(
        mhs: &MagneticHallSymbol,
        construct_type: &ConstructType,
    ) -> Operations {
        let db_prim_mag_operations = mhs.primitive_generators();
        match construct_type {
            ConstructType::Type1 | ConstructType::Type2 | ConstructType::Type3 => {
                // Reference space group: FSG
                // -> Remove time-reversal parts
                db_prim_mag_operations
                    .iter()
                    .map(|mops| mops.operation.clone())
                    .collect()
            }
            ConstructType::Type4 => {
                // Reference space group: XSG
                // -> Remove anti-translation operation
                let identity = Rotation::identity();
                db_prim_mag_operations
                    .iter()
                    .filter_map(|mops| {
                        // Here we assume the magnetic Hall symbol composes of XSG generators and one anti-translation operation
                        if mops.time_reversal {
                            assert_eq!(mops.operation.rotation, identity);
                            None
                        } else {
                            Some(mops.operation.clone())
                        }
                    })
                    .collect()
            }
        }
    }

    /// Search for origin_shift such that (trans_mat, origin_shift) transforms `prim_mag_operations` into <db_prim_mag_generators>
    /// TODO: unify with identify/space_group.rs::match_origin_shift
    fn match_origin_shift(
        prim_mag_operations: &MagneticOperations,
        trans_mat: &UnimodularLinear,
        db_prim_mag_generators: &MagneticOperations,
        epsilon: f64,
    ) -> Option<Translation> {
        let new_prim_mag_operations = UnimodularTransformation::from_linear(*trans_mat)
            .transform_magnetic_operations(prim_mag_operations);
        let mut hm_translations = HashMap::new();
        for mops in new_prim_mag_operations.iter() {
            hm_translations.insert(
                (mops.operation.rotation, mops.time_reversal),
                mops.operation.translation,
            );
        }

        let mut a = OMatrix::<i32, Dyn, U3>::zeros(3 * db_prim_mag_generators.len());
        let mut b = OVector::<f64, Dyn>::zeros(3 * db_prim_mag_generators.len());
        for (k, mops) in db_prim_mag_generators.iter().enumerate() {
            let target_translation =
                hm_translations.get(&(mops.operation.rotation, mops.time_reversal))?;

            let ak = mops.operation.rotation - Rotation::identity();
            let bk = mops.operation.translation - target_translation;
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
}

fn identify_reference_space_group(
    prim_mag_operations: &MagneticOperations,
    epsilon: f64,
) -> Option<(Operations, ConstructType)> {
    let prim_xsg = primitive_maximal_space_subgroup_from_magnetic_space_group(prim_mag_operations);
    let (fsg, is_type2) =
        family_space_group_from_magnetic_space_group(prim_mag_operations, epsilon);

    if (prim_mag_operations.len() % prim_xsg.len() != 0)
        || (prim_mag_operations.len() % fsg.len() != 0)
    {
        debug!("Input magnetic operations are incomplete.");
        return None;
    }

    let construct_type = match (prim_mag_operations.len() / prim_xsg.len(), is_type2) {
        (1, false) => ConstructType::Type1,
        (2, true) => ConstructType::Type2,
        (2, false) => {
            // Find coset representatives of MSG/XSG
            let identity = Rotation::identity();
            if prim_mag_operations
                .iter()
                .any(|mops| mops.time_reversal && mops.operation.rotation == identity)
            {
                // Anti-translation
                ConstructType::Type4
            } else {
                ConstructType::Type3
            }
        }
        _ => {
            debug!(
                "Unreachable combination: |MSG/XSG|={}, |FSG/XSG|={}",
                prim_mag_operations.len() / prim_xsg.len(),
                fsg.len() / prim_xsg.len(),
            );
            return None;
        }
    };

    // BNS setting
    let ref_spg = if construct_type == ConstructType::Type4 {
        prim_xsg
    } else {
        // For type I, II, III, `fsg` is in primitive
        fsg
    };
    Some((ref_spg, construct_type))
}

/// XSG: take only operations without time-reversal
fn primitive_maximal_space_subgroup_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
) -> Operations {
    let mut xsg = vec![];

    for mops in prim_mag_operations {
        if mops.time_reversal {
            continue;
        }
        xsg.push(mops.operation.clone());
    }
    xsg
}

/// FSG: take all operations ignoring time-reversal parts
/// Returned operations may contain duplicated rotation parts (for type-IV).
fn family_space_group_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
    epsilon: f64,
) -> (Operations, bool) {
    let mut fsg = vec![];
    let mut hm_translation = HashMap::new();
    let mut is_type2 = false;

    for mops in prim_mag_operations {
        if let Some(&other_translation) = hm_translation.get(&mops.operation.rotation) {
            let diff: Translation = mops.operation.translation - other_translation;
            if diff.iter().all(|e| (e - e.round()).abs() < epsilon) {
                is_type2 = true;
                continue;
            }
        }

        fsg.push(mops.operation.clone());
        hm_translation.insert(mops.operation.rotation.clone(), mops.operation.translation);
    }
    (fsg, is_type2)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use test_log::test as test_with_log;

    use super::*;
    use crate::base::Transformation;
    use crate::data::{MagneticHallSymbol, NUM_MAGNETIC_SPACE_GROUP_TYPES};

    fn get_prim_mag_operations(uni_number: UNINumber) -> MagneticOperations {
        let mhs = MagneticHallSymbol::from_uni_number(uni_number).unwrap();
        let magnetic_operations = mhs.traverse();

        // conventional -> primitive
        let prim_mag_operations = Transformation::from_linear(mhs.centering.linear())
            .inverse_transform_magnetic_operations(&magnetic_operations);
        prim_mag_operations
    }

    #[rstest]
    #[case(2, ConstructType::Type2, 2, 1, 1)]
    #[case(1594, ConstructType::Type1, 48, 48, 48)]
    #[case(1595, ConstructType::Type2, 96, 48, 48)]
    #[case(1596, ConstructType::Type3, 48, 24, 48)]
    #[case(1599, ConstructType::Type4, 96, 48, 96)] // -P 4 2 3 1abc' (UNI No. 1599)
    fn test_xsg_and_fsg(
        #[case] uni_number: UNINumber,
        #[case] construct_type: ConstructType,
        #[case] order_msg: usize,
        #[case] order_xsg: usize,
        #[case] order_fsg: usize,
    ) {
        let prim_mag_operations = get_prim_mag_operations(uni_number);
        assert_eq!(prim_mag_operations.len(), order_msg);

        let xsg = primitive_maximal_space_subgroup_from_magnetic_space_group(&prim_mag_operations);
        assert_eq!(xsg.len(), order_xsg);

        let epsilon = 1e-8;
        let (fsg, _) = family_space_group_from_magnetic_space_group(&prim_mag_operations, epsilon);
        assert_eq!(fsg.len(), order_fsg);

        let (_, construct_type_actual) =
            identify_reference_space_group(&prim_mag_operations, epsilon).unwrap();
        assert_eq!(construct_type_actual, construct_type);
    }

    #[test_with_log]
    fn test_identify_magnetic_space_group() {
        // for uni_number in 1..=NUM_MAGNETIC_SPACE_GROUP_TYPES {
        for uni_number in 3..=NUM_MAGNETIC_SPACE_GROUP_TYPES {
            dbg!(uni_number);
            let prim_mag_operations = get_prim_mag_operations(uni_number as UNINumber);
            let magnetic_space_group = MagneticSpaceGroup::new(&prim_mag_operations, 1e-8).unwrap();
            assert_eq!(magnetic_space_group.uni_number, uni_number as UNINumber);
        }
    }
}
