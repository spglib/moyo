use crate::base::{MagneticOperation, MagneticOperations, MoyoError, Operation, Operations};
use crate::data::hall_symbol::{HallSymbol, LayerHallSymbol, MagneticHallSymbol};
use crate::data::hall_symbol_database::{Number, hall_symbol_entry};
use crate::data::layer::{LayerNumber, LayerSetting};
use crate::data::magnetic::{UNINumber, magnetic_hall_symbol_entry};
use crate::data::setting::Setting;

/// Return symmetry operations for the given space-group ITA `number`.
///
/// `setting` selects the Hall-number convention. Pass `Setting::Spglib` for the spglib default.
/// If `primitive` is true, operations are returned for the primitive cell; otherwise for the
/// conventional cell.
pub fn operations_from_number(
    number: Number,
    setting: Setting,
    primitive: bool,
) -> Result<Operations, MoyoError> {
    let hall_number = match setting {
        Setting::HallNumber(hall_number) => hall_number,
        _ => *setting
            .hall_numbers()
            .get((number - 1) as usize)
            .ok_or(MoyoError::UnknownNumberError)?,
    };
    let entry = hall_symbol_entry(hall_number).ok_or(MoyoError::UnknownHallNumberError)?;
    let hs = HallSymbol::new(entry.hall_symbol).ok_or(MoyoError::HallSymbolParsingError)?;

    let mut operations = vec![];
    if primitive {
        for operation in hs.primitive_traverse().into_iter() {
            operations.push(operation)
        }
    } else {
        let coset = hs.traverse();

        let lattice_points = hs.centering.lattice_points();
        for t1 in lattice_points.iter() {
            for operation2 in coset.iter() {
                // (E, t1) (r2, t2) = (r2, t1 + t2)
                let t12 = (t1 + operation2.translation).map(|e| e % 1.);
                operations.push(Operation::new(operation2.rotation, t12));
            }
        }
    }

    Ok(operations)
}

/// Return symmetry operations for the given layer-group `number`.
///
/// `setting` selects the Hall-number convention. Pass `LayerSetting::Spglib` for the spglib default.
/// If `primitive` is true, operations are returned for the primitive cell; otherwise for the
/// conventional cell.
pub fn operations_from_layer_number(
    number: LayerNumber,
    setting: LayerSetting,
    primitive: bool,
) -> Result<Operations, MoyoError> {
    let hall_number = match setting {
        LayerSetting::HallNumber(hall_number) => hall_number,
        _ => *setting
            .hall_numbers()
            .get((number - 1) as usize)
            .ok_or(MoyoError::UnknownNumberError)?,
    };
    let lhs =
        LayerHallSymbol::from_hall_number(hall_number).ok_or(MoyoError::UnknownHallNumberError)?;

    let mut operations = vec![];
    if primitive {
        for operation in lhs.primitive_traverse().into_iter() {
            operations.push(operation)
        }
    } else {
        let coset = lhs.traverse();

        let lattice_points = lhs.centering().lattice_points();
        for t1 in lattice_points.iter() {
            for operation2 in coset.iter() {
                // (E, t1) (r2, t2) = (r2, t1 + t2)
                let t12 = (t1 + operation2.translation).map(|e| e.rem_euclid(1.0));
                operations.push(Operation::new(operation2.rotation, t12));
            }
        }
    }

    Ok(operations)
}

/// Return magnetic symmetry operations for the given UNI (and BNS) `uni_number`.
///
/// If `primitive` is true, operations are returned for the primitive cell; otherwise for the
/// conventional cell.
pub fn magnetic_operations_from_uni_number(
    uni_number: UNINumber,
    primitive: bool,
) -> Result<MagneticOperations, MoyoError> {
    let entry = magnetic_hall_symbol_entry(uni_number).ok_or(MoyoError::UnknownUNINumberError)?;
    let mhs = MagneticHallSymbol::new(entry.magnetic_hall_symbol)
        .ok_or(MoyoError::HallSymbolParsingError)?;

    let mut magnetic_operations = vec![];
    if primitive {
        for mops in mhs.primitive_traverse().into_iter() {
            magnetic_operations.push(mops)
        }
    } else {
        let coset = mhs.traverse();

        let lattice_points = mhs.centering.lattice_points();
        for t1 in lattice_points.iter() {
            for mops2 in coset.iter() {
                // (E, t1) (r2, t2) = (r2, t1 + t2)
                let t12 = (t1 + mops2.operation.translation).map(|e| e % 1.);
                magnetic_operations.push(MagneticOperation::new(
                    mops2.operation.rotation,
                    t12,
                    mops2.time_reversal,
                ));
            }
        }
    }

    Ok(magnetic_operations)
}

#[cfg(test)]
mod tests {
    use nalgebra::vector;

    use super::*;
    use crate::base::Position;

    fn unique_sites(position: &Position, operations: &Operations) -> Vec<Position> {
        let mut sites: Vec<Position> = vec![];
        for operation in operations.iter() {
            let new_site = operation.rotation.map(|e| e as f64) * position + operation.translation;
            let mut overlap = false;
            for site in sites.iter() {
                let mut diff = site - new_site;
                diff -= diff.map(|x| x.round());
                if diff.iter().all(|x| x.abs() < 1e-4) {
                    overlap = true;
                    break;
                }
            }
            if !overlap {
                sites.push(new_site);
            }
        }
        sites
    }

    #[test]
    fn test_operations_from_number() {
        {
            // C2/c
            let operations = operations_from_number(15, Setting::Spglib, false).unwrap();
            let x = 0.1234;
            let y = 0.5678;
            let z = 0.9012;
            assert!(unique_sites(&vector![0.0, 0.0, 0.0], &operations).len() == 4);
            assert!(unique_sites(&vector![0.0, 0.5, 0.0], &operations).len() == 4);
            assert!(unique_sites(&vector![0.25, 0.25, 0.0], &operations).len() == 4);
            assert!(unique_sites(&vector![0.25, 0.25, 0.5], &operations).len() == 4);
            assert!(unique_sites(&vector![0.0, y, 0.25], &operations).len() == 4);
            assert!(unique_sites(&vector![x, y, z], &operations).len() == 8);
        }
        {
            // Pm-3m
            let operations = operations_from_number(221, Setting::Spglib, false).unwrap();
            let x = 0.1234;
            let y = 0.5678;
            let z = 0.9012;
            assert!(operations.len() == 48);
            assert!(unique_sites(&vector![0.0, 0.0, 0.0], &operations).len() == 1);
            assert!(unique_sites(&vector![0.5, 0.5, 0.5], &operations).len() == 1);
            assert!(unique_sites(&vector![0.0, 0.5, 0.5], &operations).len() == 3);
            assert!(unique_sites(&vector![0.5, 0.0, 0.0], &operations).len() == 3);
            assert!(unique_sites(&vector![x, 0.0, 0.0], &operations).len() == 6);
            assert!(unique_sites(&vector![x, 0.5, 0.5], &operations).len() == 6);
            assert!(unique_sites(&vector![x, x, x], &operations).len() == 8);
            assert!(unique_sites(&vector![x, 0.5, 0.0], &operations).len() == 12);
            assert!(unique_sites(&vector![0.0, y, y], &operations).len() == 12);
            assert!(unique_sites(&vector![0.5, y, y], &operations).len() == 12);
            assert!(unique_sites(&vector![0.0, y, z], &operations).len() == 24);
            assert!(unique_sites(&vector![0.5, y, z], &operations).len() == 24);
            assert!(unique_sites(&vector![x, x, z], &operations).len() == 24);
            assert!(unique_sites(&vector![x, y, z], &operations).len() == 48);
        }
        {
            // Im-3m
            let operations = operations_from_number(229, Setting::Spglib, false).unwrap();
            let x = 0.1234;
            let y = 0.5678;
            let z = 0.9012;
            assert!(unique_sites(&vector![0.0, 0.0, 0.0], &operations).len() == 2);
            assert!(unique_sites(&vector![0.0, 0.5, 0.5], &operations).len() == 6);
            assert!(unique_sites(&vector![0.25, 0.25, 0.25], &operations).len() == 8);
            assert!(unique_sites(&vector![0.25, 0.0, 0.5], &operations).len() == 12);
            assert!(unique_sites(&vector![x, 0.0, 0.0], &operations).len() == 12);
            assert!(unique_sites(&vector![x, x, x], &operations).len() == 16);
            assert!(unique_sites(&vector![x, 0.0, 0.5], &operations).len() == 24);
            assert!(unique_sites(&vector![0.0, y, y], &operations).len() == 24);
            assert!(unique_sites(&vector![0.25, y, 0.5 - y], &operations).len() == 48);
            assert!(unique_sites(&vector![0.0, y, z], &operations).len() == 48);
            assert!(unique_sites(&vector![x, x, z], &operations).len() == 48);
            assert!(unique_sites(&vector![x, y, z], &operations).len() == 96);
        }
        {
            // Ia-3d
            let operations = operations_from_number(230, Setting::Spglib, false).unwrap();
            let x = 0.1234;
            let y = 0.5678;
            let z = 0.9012;

            assert!(unique_sites(&vector![0.0, 0.0, 0.0], &operations).len() == 16);
            assert!(unique_sites(&vector![0.125, 0.125, 0.125], &operations).len() == 16);
            assert!(unique_sites(&vector![0.125, 0.0, 0.25], &operations).len() == 24);
            assert!(unique_sites(&vector![0.375, 0.0, 0.25], &operations).len() == 24);
            assert!(unique_sites(&vector![x, x, x], &operations).len() == 32);
            assert!(unique_sites(&vector![x, 0.0, 0.25], &operations).len() == 48);
            assert!(unique_sites(&vector![0.125, y, 0.25 - y], &operations).len() == 48);
            assert!(unique_sites(&vector![x, y, z], &operations).len() == 96);
        }
        {
            // primitive=true
            let prim_operations = operations_from_number(230, Setting::Spglib, true).unwrap();
            assert!(prim_operations.len() == 48);
        }
    }

    #[test]
    fn test_magnetic_operations_from_uni_number() {
        // UNI: R31'_c[R3] (1242), BNS: R_I3 (146.12)
        let magnetic_operations = magnetic_operations_from_uni_number(1242, false).unwrap();
        assert_eq!(magnetic_operations.len(), 18);
        let primitive_magnetic_operations =
            magnetic_operations_from_uni_number(1242, true).unwrap();
        assert_eq!(primitive_magnetic_operations.len(), 6);
    }
}
