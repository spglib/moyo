pub mod base;
pub mod data;
pub mod dataset;

#[::safer_ffi::cfg_headers]
#[test]
pub fn generate_headers() -> ::std::io::Result<()> {
    ::safer_ffi::headers::builder()
        .to_file("moyoc.h")?
        .generate()
}
