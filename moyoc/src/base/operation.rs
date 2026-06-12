use moyo::base::{Operation, Operations};
use moyo::utils::{to_3_slice, to_3x3_slice, to_matrix3, to_vector3};

use crate::ffi::{free_slice, leak_slice};

/// Symmetry operations.
/// All pointer fields are owned by the containing `MoyoDataset` and are freed
/// by `moyo_dataset_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoOperations {
    /// Rotation parts of the `num_operations` operations
    pub rotations: *const [[i32; 3]; 3],
    /// Translation parts of the `num_operations` operations
    pub translations: *const [f64; 3],
    /// Number of symmetry operations
    pub num_operations: i32,
}

impl From<&Operations> for MoyoOperations {
    fn from(operations: &Operations) -> Self {
        let rotations = operations
            .iter()
            .map(|x| to_3x3_slice(&x.rotation))
            .collect::<Vec<_>>();
        let translations = operations
            .iter()
            .map(|x| to_3_slice(&x.translation))
            .collect::<Vec<_>>();
        Self {
            rotations: leak_slice(rotations),
            translations: leak_slice(translations),
            num_operations: operations.len() as i32,
        }
    }
}

impl From<&MoyoOperations> for Operations {
    fn from(operations: &MoyoOperations) -> Self {
        let num_operations = operations.num_operations as usize;
        let rotations = unsafe {
            std::slice::from_raw_parts(operations.rotations, num_operations)
                .iter()
                .map(to_matrix3)
                .collect::<Vec<_>>()
        };
        let translations = unsafe {
            std::slice::from_raw_parts(operations.translations, num_operations)
                .iter()
                .map(to_vector3)
                .collect::<Vec<_>>()
        };
        rotations
            .into_iter()
            .zip(translations)
            .map(|(r, t)| Operation {
                rotation: r,
                translation: t,
            })
            .collect()
    }
}

pub(crate) unsafe fn moyo_operations_free_members(operations: &MoyoOperations) {
    unsafe {
        free_slice(operations.rotations, operations.num_operations as usize);
        free_slice(operations.translations, operations.num_operations as usize);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use nalgebra::{matrix, vector};

    use moyo::base::Operation;

    #[test]
    fn test_roundtrip_moyo_operations() {
        // P3_1
        let original = vec![
            Operation::new(
                matrix![
                    1, 0, 0;
                    0, 1, 0;
                    0, 0, 1;
                ],
                vector![0.0, 0.0, 0.0],
            ),
            Operation::new(
                matrix![
                    0, -1, 0;
                    1, -1, 0;
                    0, 0, 1;
                ],
                vector![0.0, 0.0, 1.0 / 3.0],
            ),
            Operation::new(
                matrix![
                    -1, 1, 0;
                    -1, 0, 0;
                    0, 0, 1;
                ],
                vector![0.0, 0.0, 2.0 / 3.0],
            ),
        ];
        let moyoc = MoyoOperations::from(&original);
        let reconstructed = Operations::from(&moyoc);
        assert_eq!(original.len(), reconstructed.len());
        for (op1, op2) in original.iter().zip(reconstructed.iter()) {
            assert_eq!(op1.rotation, op2.rotation);
            assert_relative_eq!(op1.translation, op2.translation);
        }
        unsafe { moyo_operations_free_members(&moyoc) };
    }
}
