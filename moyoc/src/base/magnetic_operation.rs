use moyo::base::{MagneticOperation, MagneticOperations};
use moyo::utils::{to_3_slice, to_3x3_slice, to_matrix3, to_vector3};

use crate::ffi::{free_slice, leak_slice};

/// Magnetic symmetry operations.
/// All pointer fields are owned by the containing dataset and are freed
/// together with it.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoMagneticOperations {
    /// Rotation parts of the `num_operations` magnetic operations
    pub rotations: *const [[i32; 3]; 3],
    /// Translation parts of the `num_operations` magnetic operations
    pub translations: *const [f64; 3],
    /// Time-reversal parts of the `num_operations` magnetic operations
    pub time_reversals: *const bool,
    /// Number of magnetic symmetry operations
    pub num_operations: i32,
}

impl From<&MagneticOperations> for MoyoMagneticOperations {
    fn from(magnetic_operations: &MagneticOperations) -> Self {
        let rotations = magnetic_operations
            .iter()
            .map(|x| to_3x3_slice(&x.operation.rotation))
            .collect::<Vec<_>>();
        let translations = magnetic_operations
            .iter()
            .map(|x| to_3_slice(&x.operation.translation))
            .collect::<Vec<_>>();
        let time_reversals = magnetic_operations
            .iter()
            .map(|x| x.time_reversal)
            .collect::<Vec<_>>();
        Self {
            rotations: leak_slice(rotations),
            translations: leak_slice(translations),
            time_reversals: leak_slice(time_reversals),
            num_operations: magnetic_operations.len() as i32,
        }
    }
}

impl From<&MoyoMagneticOperations> for MagneticOperations {
    fn from(magnetic_operations: &MoyoMagneticOperations) -> Self {
        let num_operations = magnetic_operations.num_operations as usize;
        unsafe {
            let rotations =
                std::slice::from_raw_parts(magnetic_operations.rotations, num_operations);
            let translations =
                std::slice::from_raw_parts(magnetic_operations.translations, num_operations);
            let time_reversals =
                std::slice::from_raw_parts(magnetic_operations.time_reversals, num_operations);
            rotations
                .iter()
                .zip(translations)
                .zip(time_reversals)
                .map(|((r, t), &tr)| MagneticOperation::new(to_matrix3(r), to_vector3(t), tr))
                .collect()
        }
    }
}

pub(crate) unsafe fn moyo_magnetic_operations_free_members(
    magnetic_operations: &MoyoMagneticOperations,
) {
    unsafe {
        let num_operations = magnetic_operations.num_operations as usize;
        free_slice(magnetic_operations.rotations, num_operations);
        free_slice(magnetic_operations.translations, num_operations);
        free_slice(magnetic_operations.time_reversals, num_operations);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use nalgebra::{matrix, vector};

    #[test]
    fn test_roundtrip_moyo_magnetic_operations() {
        // P-1' generators
        let original = vec![
            MagneticOperation::new(
                matrix![
                    1, 0, 0;
                    0, 1, 0;
                    0, 0, 1;
                ],
                vector![0.0, 0.0, 0.0],
                false,
            ),
            MagneticOperation::new(
                matrix![
                    -1, 0, 0;
                    0, -1, 0;
                    0, 0, -1;
                ],
                vector![0.0, 0.0, 0.0],
                true,
            ),
        ];
        let moyoc = MoyoMagneticOperations::from(&original);
        let reconstructed = MagneticOperations::from(&moyoc);
        assert_eq!(original.len(), reconstructed.len());
        for (op1, op2) in original.iter().zip(reconstructed.iter()) {
            assert_eq!(op1.operation.rotation, op2.operation.rotation);
            assert_relative_eq!(op1.operation.translation, op2.operation.translation);
            assert_eq!(op1.time_reversal, op2.time_reversal);
        }
        unsafe { moyo_magnetic_operations_free_members(&moyoc) };
    }
}
