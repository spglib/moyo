use moyo::data::Setting;

#[repr(C)]
pub enum MoyocSetting {
    SPGLIB,
    STANDARD,
}

impl From<MoyocSetting> for Setting {
    fn from(setting: MoyocSetting) -> Self {
        match setting {
            MoyocSetting::SPGLIB => Setting::Spglib,
            MoyocSetting::STANDARD => Setting::Standard,
        }
    }
}
