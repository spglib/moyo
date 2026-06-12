use moyo::base::Lattice;
use moyo::identify::PointGroup;
use moyo::utils::{to_3x3_slice, to_matrix3};

/// Point group identified by `moyo_point_group_new` and freed by
/// `moyo_point_group_free`.
#[derive(Debug, Clone)]
#[repr(C)]
pub struct MoyoPointGroup {
    /// Number for the arithmetic crystal class (1 - 73).
    pub arithmetic_number: i32,
    /// Transformation matrix from the input primitive basis to the standardized basis.
    pub prim_trans_mat: [[i32; 3]; 3],
}

/// Identify the point group of the given primitive rotations.
///
/// - `prim_rotations`: rotation matrices of the `num_operations` operations in
///   the primitive cell.
/// - `basis`: row-wise basis vectors of the primitive lattice, `basis[i]` is
///   the `i`th basis vector. Pass NULL to assume an identity basis.
/// - Returns NULL if the identification fails or the arguments are invalid.
///   Free the result with `moyo_point_group_free`.
///
/// # Safety
/// `prim_rotations` must point to `num_operations` elements, and `basis` must
/// point to 3 basis vectors or be NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_point_group_new(
    prim_rotations: *const [[i32; 3]; 3],
    num_operations: i32,
    basis: *const [f64; 3],
) -> *mut MoyoPointGroup {
    if prim_rotations.is_null() || num_operations <= 0 {
        return std::ptr::null_mut();
    }
    let prim_rotations = unsafe {
        std::slice::from_raw_parts(prim_rotations, num_operations as usize)
            .iter()
            .map(to_matrix3)
            .collect::<Vec<_>>()
    };

    let point_group = if basis.is_null() {
        PointGroup::new(&prim_rotations)
    } else {
        let basis = unsafe { *(basis as *const [[f64; 3]; 3]) };
        PointGroup::from_lattice(&Lattice::from_basis(basis), &prim_rotations)
    };
    match point_group {
        Ok(point_group) => Box::into_raw(Box::new(MoyoPointGroup {
            arithmetic_number: point_group.arithmetic_number,
            prim_trans_mat: to_3x3_slice(&point_group.prim_trans_mat),
        })),
        Err(_) => std::ptr::null_mut(),
    }
}

/// Free a point group created by `moyo_point_group_new`. Passing NULL is a no-op.
///
/// # Safety
/// `point_group` must be a pointer returned by `moyo_point_group_new` that has
/// not been freed yet, or NULL.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn moyo_point_group_free(point_group: *mut MoyoPointGroup) {
    if point_group.is_null() {
        return;
    }
    unsafe {
        drop(Box::from_raw(point_group));
    }
}
