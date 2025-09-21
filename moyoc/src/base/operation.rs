use moyo::base::{Operation, Operations};
use moyo::utils::{to_3_slice, to_3x3_slice, to_matrix3, to_vector3};

#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoOperations {
    pub rotations: *const [[i32; 3]; 3],
    pub translations: *const [f64; 3],
    pub num_operations: usize,
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
            rotations: rotations.leak().as_ptr(),
            translations: translations.leak().as_ptr(),
            num_operations: operations.len(),
        }
    }
}

impl From<&MoyoOperations> for Operations {
    fn from(operations: &MoyoOperations) -> Self {
        let rotations = unsafe {
            std::slice::from_raw_parts(operations.rotations, operations.num_operations)
                .iter()
                .map(to_matrix3)
                .collect::<Vec<_>>()
        };
        let translations = unsafe {
            std::slice::from_raw_parts(operations.translations, operations.num_operations)
                .iter()
                .map(to_vector3)
                .collect::<Vec<_>>()
        };
        rotations
            .into_iter()
            .zip(translations.into_iter())
            .map(|(r, t)| Operation {
                rotation: r,
                translation: t,
            })
            .collect()
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn free_moyo_operations(operations: MoyoOperations) {
    unsafe {
        let _ = Vec::from_raw_parts(
            operations.rotations as *mut [[i32; 3]; 3],
            operations.num_operations as usize,
            operations.num_operations as usize,
        );
        let _ = Vec::from_raw_parts(
            operations.translations as *mut [f64; 3],
            operations.num_operations as usize,
            operations.num_operations as usize,
        );
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
        free_moyo_operations(moyoc);
    }
}
