use std::collections::HashMap;

use log::debug;

use super::normalizer::integral_normalizer;
use super::point_group::{iter_trans_mat_basis, iter_unimodular_trans_mat};
use super::rotation_type::identify_rotation_type;
use super::space_group::{match_origin_shift, SpaceGroup};
use crate::base::{
    project_rotations, MagneticOperations, MoyoError, Operation, Operations, Rotation, Translation,
    UnimodularTransformation,
};
use crate::data::{
    get_magnetic_space_group_type, hall_symbol_entry, magnetic_hall_symbol_entry, uni_number_range,
    ConstructType, HallSymbol, MagneticHallSymbol, MagneticHallSymbolEntry, Setting, UNINumber,
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

            let entry = magnetic_hall_symbol_entry(uni_number).unwrap();
            let mhs = MagneticHallSymbol::new(&entry.magnetic_hall_symbol)
                .ok_or(MoyoError::MagneticSpaceGroupTypeIdentificationError)?;
            let db_prim_mag_operations = mhs.primitive_traverse();
            let (db_ref_prim_operations, db_ref_prim_generators) =
                db_reference_space_group_primitive(&entry);

            match construct_type {
                ConstructType::Type3 => {
                    // Centralizer of FSG also keep XSG invariant. Thus, we only need to consider the normalizer of FSG up to the centralizer.
                    // TODO: precompute the normalizer
                    let db_fsg_prim_normalizer = integral_normalizer(
                        &db_ref_prim_operations,
                        &db_ref_prim_generators,
                        epsilon,
                    );
                    for corr_trans in db_fsg_prim_normalizer {
                        // new_transformation * corr_trans: primitive input -> primitive DB
                        let new_transformation = std_ref_spg.transformation.clone() * corr_trans;
                        let new_prim_mag_operations =
                            new_transformation.transform_magnetic_operations(prim_mag_operations);

                        if Self::match_prim_mag_operations(
                            &new_prim_mag_operations,
                            &db_prim_mag_operations,
                            epsilon,
                        ) {
                            debug!("Matched with UNI number {}", uni_number);
                            return Ok(Self {
                                uni_number,
                                transformation: new_transformation,
                            });
                        }
                    }
                }
                ConstructType::Type4 => {
                    // Find conjugator which transforms the anti-translation while keeping XSGs
                    let identity = Rotation::identity();
                    let original_anti_translation = prim_mag_operations
                        .iter()
                        .filter_map(|mops| {
                            if (mops.operation.rotation == identity) && mops.time_reversal {
                                Some(mops)
                            } else {
                                None
                            }
                        })
                        .nth(0)
                        .unwrap();
                    let src_translation = std_ref_spg
                        .transformation
                        .transform_magnetic_operation(original_anti_translation)
                        .operation
                        .translation;
                    // TODO: refactor filtering anti-translation
                    let dst_translation = db_prim_mag_operations
                        .iter()
                        .filter_map(|mops| {
                            if (mops.operation.rotation == identity) && mops.time_reversal {
                                Some(mops)
                            } else {
                                None
                            }
                        })
                        .nth(0)
                        .unwrap()
                        .operation
                        .translation;
                    if let Some(corr_trans) = find_conjugator_type4(
                        &db_ref_prim_generators,
                        &db_ref_prim_operations,
                        &src_translation,
                        &dst_translation,
                        epsilon,
                    ) {
                        let new_transformation = std_ref_spg.transformation.clone() * corr_trans;
                        let new_prim_mag_operations =
                            new_transformation.transform_magnetic_operations(&prim_mag_operations);
                        if Self::match_prim_mag_operations(
                            &new_prim_mag_operations,
                            &db_prim_mag_operations,
                            epsilon,
                        ) {
                            return Ok(Self {
                                uni_number,
                                transformation: new_transformation,
                            });
                        }
                    }
                }
                _ => unreachable!(),
            }
        }
        Err(MoyoError::MagneticSpaceGroupTypeIdentificationError)
    }

    pub fn reference_space_group(&self) -> SpaceGroup {
        let ref_hall_number = magnetic_hall_symbol_entry(self.uni_number)
            .unwrap()
            .reference_hall_number();
        SpaceGroup::from_hall_number_and_transformation(
            ref_hall_number,
            self.transformation.clone(),
        )
        .unwrap()
    }

    fn match_prim_mag_operations(
        prim_mag_operations1: &MagneticOperations,
        prim_mag_operations2: &MagneticOperations,
        epsilon: f64,
    ) -> bool {
        if prim_mag_operations1.len() != prim_mag_operations2.len() {
            return false;
        }

        let mut hm_translation = HashMap::new();
        for mops1 in prim_mag_operations1 {
            hm_translation.insert(
                (mops1.operation.rotation.clone(), mops1.time_reversal),
                mops1.operation.translation,
            );
        }

        for mops2 in prim_mag_operations2 {
            if let Some(translation1) =
                hm_translation.get(&(mops2.operation.rotation.clone(), mops2.time_reversal))
            {
                let diff = mops2.operation.translation - translation1;
                if !diff.iter().all(|e| (e - e.round()).abs() < epsilon) {
                    return false;
                }
            } else {
                return false;
            }
        }
        true
    }
}

fn identify_reference_space_group(
    prim_mag_operations: &MagneticOperations,
    epsilon: f64,
) -> Option<(Operations, ConstructType)> {
    let (prim_xsg, _) =
        primitive_maximal_space_subgroup_from_magnetic_space_group(prim_mag_operations);
    let (fsg, is_type2, _) =
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
pub fn primitive_maximal_space_subgroup_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
) -> (Operations, Vec<bool>) {
    let mut xsg = vec![];
    let mut contained = vec![false; prim_mag_operations.len()];

    for (i, mops) in prim_mag_operations.iter().enumerate() {
        if mops.time_reversal {
            continue;
        }
        xsg.push(mops.operation.clone());
        contained[i] = true;
    }
    (xsg, contained)
}

/// FSG: take all operations ignoring time-reversal parts
/// Returned operations may contain duplicated rotation parts (for type-IV).
pub fn family_space_group_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
    epsilon: f64,
) -> (Operations, bool, Vec<bool>) {
    let mut fsg = vec![];
    let mut hm_translation = HashMap::new();
    let mut contained = vec![false; prim_mag_operations.len()];
    let mut is_type2 = false;

    for (i, mops) in prim_mag_operations.iter().enumerate() {
        if let Some(&other_translation) = hm_translation.get(&mops.operation.rotation) {
            let diff: Translation = mops.operation.translation - other_translation;
            if diff.iter().all(|e| (e - e.round()).abs() < epsilon) {
                is_type2 = true;
                continue;
            }
        }

        fsg.push(mops.operation.clone());
        hm_translation.insert(mops.operation.rotation.clone(), mops.operation.translation);
        contained[i] = true;
    }
    (fsg, is_type2, contained)
}

/// Return operations and generators of the reference space group of magnetic space group `entry`.
/// This function assumes the magnetic Hall symbol is extended from the Hall symbol in the standard setting.
fn db_reference_space_group_primitive(entry: &MagneticHallSymbolEntry) -> (Operations, Operations) {
    let ref_hall_entry = hall_symbol_entry(entry.reference_hall_number()).unwrap();
    let ref_hall_symbol = HallSymbol::new(&ref_hall_entry.hall_symbol).unwrap();
    let ref_prim_operations = ref_hall_symbol.primitive_traverse();
    let identity = Rotation::identity();

    let mut ref_prim_generators = ref_hall_symbol
        .primitive_generators()
        .into_iter()
        .filter(|ops| ops.rotation != identity) // In primitive, if rotation part is identity, it is a pure translation
        .collect::<Vec<_>>();
    if ref_prim_generators.is_empty() {
        ref_prim_generators.push(Operation::identity());
    }

    (ref_prim_operations, ref_prim_generators)
}

/// Find a unimodular transformation that transforms (E, src_translation) to (E, dst_translation) while keeping `stabilized_prim_operations`, which are generated by `stabilized_prim_generators`.
fn find_conjugator_type4(
    stabilized_prim_generators: &Operations,
    stabilized_prim_operations: &Operations,
    src_translation: &Translation,
    dst_translation: &Translation,
    epsilon: f64,
) -> Option<UnimodularTransformation> {
    let stabilized_prim_rotations = project_rotations(stabilized_prim_operations);
    let stabilized_prim_rotation_generators = project_rotations(stabilized_prim_generators);

    let rotation_types = stabilized_prim_rotations
        .iter()
        .map(identify_rotation_type)
        .collect::<Vec<_>>();
    for trans_mat_basis in iter_trans_mat_basis(
        stabilized_prim_rotations,
        rotation_types,
        stabilized_prim_rotation_generators,
    ) {
        for prim_trans_mat in iter_unimodular_trans_mat(trans_mat_basis) {
            // (P, p)^-1 (E, c_src) (P, p) = (P^-1, -P^-1 p) (P, p + c_src) = (E, P^-1 c_src) == (E, c_dst)
            let diff = prim_trans_mat.map(|e| e as f64) * dst_translation - src_translation;
            if !diff.iter().all(|e| (e - e.round()).abs() < epsilon) {
                continue;
            }

            if let Some(origin_shift) = match_origin_shift(
                stabilized_prim_operations,
                &prim_trans_mat,
                stabilized_prim_generators,
                epsilon,
            ) {
                return Some(UnimodularTransformation::new(prim_trans_mat, origin_shift));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use test_log::test as test_with_log;

    use super::*;
    use crate::data::{
        magnetic_hall_symbol_entry, MagneticHallSymbol, NUM_MAGNETIC_SPACE_GROUP_TYPES,
    };

    fn get_prim_mag_operations(uni_number: UNINumber) -> MagneticOperations {
        let mhs = MagneticHallSymbol::from_uni_number(uni_number).unwrap();
        mhs.primitive_traverse()
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

        let (xsg, _) =
            primitive_maximal_space_subgroup_from_magnetic_space_group(&prim_mag_operations);
        assert_eq!(xsg.len(), order_xsg);

        let epsilon = 1e-8;
        let (fsg, _, _) =
            family_space_group_from_magnetic_space_group(&prim_mag_operations, epsilon);
        assert_eq!(fsg.len(), order_fsg);

        let (_, construct_type_actual) =
            identify_reference_space_group(&prim_mag_operations, epsilon).unwrap();
        assert_eq!(construct_type_actual, construct_type);
    }

    // Check generators of reference space group by two methods:
    // 1. From the magnetic Hall symbol
    // 2. From the Hall symbol with the corresponding Hall number
    #[test_with_log]
    fn test_db_reference_space_group_primitive_generators() {
        for uni_number in 1..=NUM_MAGNETIC_SPACE_GROUP_TYPES {
            let entry = magnetic_hall_symbol_entry(uni_number as UNINumber).unwrap();
            let (_, actual) = db_reference_space_group_primitive(&entry);

            let mhs = MagneticHallSymbol::new(&entry.magnetic_hall_symbol).unwrap();
            let identity = Rotation::identity();
            let mut expect: Operations = match entry.construct_type() {
                ConstructType::Type1 | ConstructType::Type2 => mhs
                    .primitive_generators()
                    .iter()
                    .filter_map(|mops| {
                        if mops.time_reversal || mops.operation.rotation == identity {
                            // Ignore 1' for Type2
                            None
                        } else {
                            Some(mops.operation.clone())
                        }
                    })
                    .collect(),
                ConstructType::Type3 => mhs
                    .primitive_generators()
                    .iter()
                    .map(|mops| mops.operation.clone()) // Ignore time-reversal parts
                    .collect(),
                ConstructType::Type4 => mhs
                    .primitive_generators()
                    .iter()
                    .filter_map(|mops| {
                        if mops.operation.rotation == identity {
                            // Ignore anti-translation
                            None
                        } else {
                            Some(mops.operation.clone())
                        }
                    })
                    .collect(),
            };
            if expect.is_empty() {
                expect.push(Operation::identity());
            }

            assert_eq!(actual.len(), expect.len());
            let mut hm_translation = HashMap::new();
            for ops1 in actual.iter() {
                hm_translation.insert(ops1.rotation.clone(), ops1.translation);
            }
            for ops2 in expect.iter() {
                let translation1 = hm_translation.get(&ops2.rotation).unwrap();
                let diff = ops2.translation - translation1;
                assert_relative_eq!(diff.map(|e| (e - e.round().abs())).max(), 0.0);
            }
        }
    }

    #[test_with_log]
    fn test_identify_magnetic_space_group() {
        for uni_number in 1..=NUM_MAGNETIC_SPACE_GROUP_TYPES {
            let prim_mag_operations = get_prim_mag_operations(uni_number as UNINumber);
            let magnetic_space_group = MagneticSpaceGroup::new(&prim_mag_operations, 1e-8).unwrap();
            assert_eq!(magnetic_space_group.uni_number, uni_number as UNINumber);
        }
    }
}
