use moyo::base::Operation;
use safer_ffi;
use safer_ffi::prelude::*;

#[derive_ReprC]
#[repr(C)]
#[derive(Debug)]
pub struct MoyocOperation {
    pub rotation: [[i32; 3]; 3],
    pub translation: [f64; 3],
}

impl From<Operation> for MoyocOperation {
    fn from(operations: Operation) -> Self {
        let r = operations.rotation;
        let rotation = [
            [r[(0, 0)], r[(0, 1)], r[(0, 2)]],
            [r[(1, 0)], r[(1, 1)], r[(1, 2)]],
            [r[(2, 0)], r[(2, 1)], r[(2, 2)]],
        ];
        let t = operations.translation;
        let translation = [t[0], t[1], t[2]];
        MoyocOperation {
            rotation,
            translation,
        }
    }
}
