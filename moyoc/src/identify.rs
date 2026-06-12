mod layer_group;
mod magnetic_space_group;
mod point_group;
mod space_group;

pub use layer_group::{MoyoLayerGroup, moyo_layer_group_free, moyo_layer_group_new};
pub use magnetic_space_group::{
    MoyoMagneticSpaceGroup, moyo_magnetic_space_group_free, moyo_magnetic_space_group_new,
};
pub use point_group::{MoyoPointGroup, moyo_point_group_free, moyo_point_group_new};
pub use space_group::{MoyoSpaceGroup, moyo_space_group_free, moyo_space_group_new};

/// Default numerical tolerance for matching translations, following moyopy.
const DEFAULT_EPSILON: f64 = 1e-4;

pub(crate) fn epsilon_or_default(epsilon: f64) -> f64 {
    if epsilon < 0.0 {
        DEFAULT_EPSILON
    } else {
        epsilon
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{
        MoyoLayerSetting, MoyoSetting, moyo_magnetic_operations_free,
        moyo_magnetic_operations_from_uni_number, moyo_operations_free,
        moyo_operations_from_layer_number, moyo_operations_from_number,
    };

    #[test]
    fn test_moyo_point_group_new() {
        // P3c1 -> arithmetic crystal class 45 (3m1P)
        let operations = moyo_operations_from_number(158, MoyoSetting::Standard, 0, true);
        assert!(!operations.is_null());
        unsafe {
            let point_group = moyo_point_group_new(
                (*operations).rotations,
                (*operations).num_operations,
                std::ptr::null(),
            );
            assert!(!point_group.is_null());
            assert_eq!((*point_group).arithmetic_number, 45);
            moyo_point_group_free(point_group);
            moyo_operations_free(operations);
        }

        // Invalid input
        let null = unsafe { moyo_point_group_new(std::ptr::null(), 0, std::ptr::null()) };
        assert!(null.is_null());
    }

    #[test]
    fn test_moyo_space_group_new() {
        // P31c (No. 159)
        let operations = moyo_operations_from_number(159, MoyoSetting::Standard, 0, true);
        assert!(!operations.is_null());
        unsafe {
            let space_group = moyo_space_group_new(
                (*operations).rotations,
                (*operations).translations,
                (*operations).num_operations,
                std::ptr::null(),
                MoyoSetting::Standard,
                0,
                -1.0,
            );
            assert!(!space_group.is_null());
            assert_eq!((*space_group).number, 159);
            moyo_space_group_free(space_group);
            moyo_operations_free(operations);
        }
    }

    #[test]
    fn test_moyo_layer_group_new() {
        // p6 (LG 65), with and without an explicit identity basis
        let operations = moyo_operations_from_layer_number(65, MoyoLayerSetting::Standard, 0, true);
        assert!(!operations.is_null());
        let identity_basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        unsafe {
            for basis in [std::ptr::null(), identity_basis.as_ptr()] {
                let layer_group = moyo_layer_group_new(
                    (*operations).rotations,
                    (*operations).translations,
                    (*operations).num_operations,
                    basis,
                    MoyoLayerSetting::Standard,
                    0,
                    -1.0,
                );
                assert!(!layer_group.is_null());
                assert_eq!((*layer_group).number, 65);
                assert!((1..=116).contains(&(*layer_group).hall_number));
                moyo_layer_group_free(layer_group);
            }
            moyo_operations_free(operations);
        }
    }

    #[test]
    fn test_moyo_magnetic_space_group_new() {
        // R31'_c (UNI 1242)
        let magnetic_operations = moyo_magnetic_operations_from_uni_number(1242, true);
        assert!(!magnetic_operations.is_null());
        unsafe {
            let magnetic_space_group = moyo_magnetic_space_group_new(
                (*magnetic_operations).rotations,
                (*magnetic_operations).translations,
                (*magnetic_operations).time_reversals,
                (*magnetic_operations).num_operations,
                std::ptr::null(),
                -1.0,
            );
            assert!(!magnetic_space_group.is_null());
            assert_eq!((*magnetic_space_group).uni_number, 1242);
            moyo_magnetic_space_group_free(magnetic_space_group);
            moyo_magnetic_operations_free(magnetic_operations);
        }
    }
}
