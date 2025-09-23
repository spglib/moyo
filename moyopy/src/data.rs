mod arithmetic_crystal_class;
mod centering;
mod hall_symbol;
mod magnetic_space_group_type;
mod setting;
mod space_group_type;

pub use arithmetic_crystal_class::PyArithmeticCrystalClass;
pub use centering::PyCentering;
pub use hall_symbol::PyHallSymbolEntry;
pub use magnetic_space_group_type::PyMagneticSpaceGroupType;
pub use setting::PySetting;
pub use space_group_type::PySpaceGroupType;

use pyo3::prelude::*;

use super::base::{PyMagneticOperations, PyMoyoError, PyOperations};
use moyo::base::{MagneticOperation, MoyoError, Operation};
use moyo::data::{
    HallSymbol, MagneticHallSymbol, Number, Setting, UNINumber, hall_symbol_entry,
    magnetic_hall_symbol_entry,
};

#[pyfunction]
#[pyo3(signature = (number, *, setting=None, primitive=false))]
pub fn operations_from_number(
    number: Number,
    setting: Option<PySetting>,
    primitive: bool,
) -> Result<PyOperations, PyMoyoError> {
    let setting = if let Some(setting) = setting {
        setting
    } else {
        PySetting(Setting::Spglib)
    };
    let hall_number = match setting.0 {
        Setting::HallNumber(hall_number) => hall_number,
        Setting::Spglib | Setting::Standard => *setting
            .0
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

    Ok(PyOperations::from(operations))
}

#[pyfunction]
#[pyo3(signature = (uni_number, *, primitive=false))]
pub fn magnetic_operations_from_uni_number(
    uni_number: UNINumber,
    primitive: bool,
) -> Result<PyMagneticOperations, PyMoyoError> {
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

    Ok(PyMagneticOperations::from(magnetic_operations))
}

#[cfg(test)]
mod tests {
    use nalgebra::vector;

    use super::*;
    use moyo::base::{Operations, Position};

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
            let operations = operations_from_number(15, None, false).unwrap();
            let operations = Operations::from(operations);
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
            let operations = operations_from_number(221, None, false).unwrap();
            let operations = Operations::from(operations);
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
            let operations = operations_from_number(229, None, false).unwrap();
            let operations = Operations::from(operations);
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
            let operations = operations_from_number(230, None, false).unwrap();
            let operations = Operations::from(operations);
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
            let prim_operations = operations_from_number(230, None, true).unwrap();
            assert!(prim_operations.num_operations() == 48);
        }
    }

    #[test]
    fn test_magnetic_operations_from_uni_number() {
        {
            // UNI: R31'_c[R3] (1242), BNS: R_I3 (146.12)
            let magnetic_operations = magnetic_operations_from_uni_number(1242, false).unwrap();
            assert_eq!(magnetic_operations.num_operations(), 18);
            let primitive_magnetic_operations =
                magnetic_operations_from_uni_number(1242, true).unwrap();
            assert_eq!(primitive_magnetic_operations.num_operations(), 6);
        }
    }
}
