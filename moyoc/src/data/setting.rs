/// Preference for the setting of the space group.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum MoyoSetting {
    /// Hall number specified by the `hall_number` argument
    HallNumber,
    /// The setting of the smallest Hall number
    Spglib,
    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral,
    /// and origin choice 2 for centrosymmetric space groups
    Standard,
}
