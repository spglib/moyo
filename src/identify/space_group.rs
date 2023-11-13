pub enum Setting {
    HallNumber(i32),
    /// The setting of the smallest Hall number
    Spglib,
    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral, and origin choice 2 for centrosymmetric space groups
    Standard,
}
