use std::collections::HashSet;

use log::debug;

use super::space_group::SpaceGroup;
use crate::base::{MagneticOperations, MoyoError, Operations, Rotation, UnimodularTransformation};
use crate::data::{uni_number_range, ConstructType, Setting, UNINumber};

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
        let prim_xsg = maximal_space_subgroup_from_magnetic_space_group(prim_mag_operations)
            .ok_or(MoyoError::MagneticSpaceGroupTypeIdentificationError)?;
        let prim_fsg = family_space_group_from_magnetic_space_group(prim_mag_operations)
            .ok_or(MoyoError::MagneticSpaceGroupTypeIdentificationError)?;

        let construct_type = match (
            prim_mag_operations.len() % prim_xsg.len(),
            prim_fsg.len() % prim_xsg.len(),
        ) {
            (1, 1) => ConstructType::Type1,
            (2, 1) => ConstructType::Type2,
            (2, 2) => {
                // Find coset representatives of MSG/XSG
                let identify = Rotation::identity();
                if prim_mag_operations
                    .iter()
                    .any(|mops| mops.time_reversal && mops.operation.rotation == identify)
                {
                    // Anti-translation
                    ConstructType::Type4
                } else {
                    ConstructType::Type3
                }
            }
            _ => {
                return Err(MoyoError::MagneticSpaceGroupTypeIdentificationError);
            }
        };

        // BNS setting
        // std_ref_spg.transformation: primitive input -> primitive BNS setting
        let ref_spg = if construct_type == ConstructType::Type4 {
            prim_xsg
        } else {
            prim_fsg
        };
        let std_ref_spg = SpaceGroup::new(&ref_spg, Setting::Standard, epsilon)?;
        let tmp_prim_mag_operations = std_ref_spg
            .transformation
            .transform_magnetic_operations(&prim_mag_operations);

        let uni_number_range = uni_number_range(std_ref_spg.number)
            .ok_or(MoyoError::SpaceGroupTypeIdentificationError)?;
        for uni_number in uni_number_range {
            todo!()
        }

        todo!()
    }
}

fn maximal_space_subgroup_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
) -> Option<Operations> {
    let mut xsg = vec![];
    let mut visited = HashSet::new();

    for mops in prim_mag_operations {
        if mops.time_reversal {
            continue;
        }

        if visited.contains(&mops.operation.rotation) {
            debug!("Input magnetic operations are not in primitive.");
            return None;
        }
        xsg.push(mops.operation.clone());
        visited.insert(mops.operation.rotation.clone());
    }

    if (xsg.len() != prim_mag_operations.len()) && (xsg.len() * 2 != prim_mag_operations.len()) {
        debug!("Input magnetic operations are incomplete.");
        return None;
    }
    Some(xsg)
}

fn family_space_group_from_magnetic_space_group(
    prim_mag_operations: &MagneticOperations,
) -> Option<Operations> {
    let mut fsg = vec![];
    let mut visited = HashSet::new();

    for mops in prim_mag_operations {
        if visited.contains(&mops.operation.rotation) {
            continue;
        }
        fsg.push(mops.operation.clone());
        visited.insert(mops.operation.rotation.clone());
    }

    if (fsg.len() != prim_mag_operations.len()) && (fsg.len() * 2 != prim_mag_operations.len()) {
        debug!("Input magnetic operations are incomplete.");
        return None;
    }
    Some(fsg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Transformation;
    use crate::data::MagneticHallSymbol;

    #[test]
    fn test_type4() {
        // -P 4 2 3 1abc' (UNI No. 1599)
        // P_I m-3m (BNS No. 221.97)
        let mhs = MagneticHallSymbol::from_uni_number(1599).unwrap();
        let magnetic_operations = mhs.traverse();

        // conventional -> primitive
        let prim_mag_operations = Transformation::from_linear(mhs.centering.linear())
            .inverse_transform_magnetic_operations(&magnetic_operations);

        let magnetic_space_group = MagneticSpaceGroup::new(&prim_mag_operations, 1e-8).unwrap();
    }
}
