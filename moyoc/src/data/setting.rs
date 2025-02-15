use ::safer_ffi::prelude::*;
use moyo::data::Setting;

#[derive_ReprC]
#[repr(u8)]
pub enum MoyocSetting {
    Spglib,
    Standard,
}

impl From<MoyocSetting> for Setting {
    fn from(setting: MoyocSetting) -> Self {
        match setting {
            MoyocSetting::Spglib => Setting::Spglib,
            MoyocSetting::Standard => Setting::Standard,
        }
    }
}
