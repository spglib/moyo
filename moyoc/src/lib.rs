#[cfg(test)]
#[macro_use]
extern crate approx;

use core::ffi::c_char;

pub mod base;
pub mod data;
pub mod dataset;

pub use base::{
    MoyoCell, MoyoCollinearMagneticCell, MoyoMagneticOperations, MoyoNonCollinearMagneticCell,
    MoyoOperations,
};
pub use data::{
    MoyoLayerSetting, MoyoSetting, moyo_magnetic_operations_free,
    moyo_magnetic_operations_from_uni_number, moyo_operations_free,
    moyo_operations_from_layer_number, moyo_operations_from_number,
};
pub use dataset::{
    MoyoCollinearMagneticDataset, MoyoDataset, MoyoLayerDataset, MoyoNonCollinearMagneticDataset,
    moyo_collinear_magnetic_dataset_free, moyo_collinear_magnetic_dataset_new, moyo_dataset_free,
    moyo_dataset_new, moyo_layer_dataset_free, moyo_layer_dataset_new,
    moyo_noncollinear_magnetic_dataset_free, moyo_noncollinear_magnetic_dataset_new,
};

/// Return the version string of moyoc (e.g. "0.11.0").
/// The returned string is statically allocated and must not be freed.
#[unsafe(no_mangle)]
pub extern "C" fn moyo_version() -> *const c_char {
    concat!(env!("CARGO_PKG_VERSION"), "\0").as_ptr() as *const c_char
}

pub(crate) mod ffi {
    use core::ffi::c_char;
    use std::ffi::CString;

    /// Leak a vector as a raw pointer to be owned by the C side.
    pub fn leak_slice<T>(v: Vec<T>) -> *const T {
        Box::leak(v.into_boxed_slice()).as_ptr()
    }

    /// Free a pointer created by `leak_slice` with the same length.
    pub unsafe fn free_slice<T>(ptr: *const T, len: usize) {
        if ptr.is_null() {
            return;
        }
        unsafe {
            drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
                ptr as *mut T,
                len,
            )));
        }
    }

    /// Leak a string as a NUL-terminated C string to be owned by the C side.
    pub fn leak_cstring(s: String) -> *const c_char {
        CString::new(s).expect("CString::new failed").into_raw()
    }

    /// Free a pointer created by `leak_cstring`.
    pub unsafe fn free_cstring(ptr: *const c_char) {
        if ptr.is_null() {
            return;
        }
        unsafe {
            drop(CString::from_raw(ptr as *mut c_char));
        }
    }
}
