/// Preference for the Hall setting of the layer group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum MoyoLayerSetting {
    /// Layer Hall number specified by the `hall_number` argument (1 - 116)
    HallNumber,
    /// The setting of the smallest layer Hall number for each layer group
    /// (origin choice 1 for the centrosymmetric layer groups 52, 62, 64)
    Spglib,
    /// BCS standard choice (de la Flor et al., Acta Cryst. A77, 559-571 (2021)):
    /// origin choice 2 for the centrosymmetric layer groups 52, 62, 64
    Standard,
}
