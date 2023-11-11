pub enum CrystalFamily {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Hexagonal,
    Cubic,
}

pub enum LatticeSystem {
    // (lattice system, holohedry)
    Triclinic,    // -1
    Monoclinic,   // 2/m
    Orthorhombic, // mmm
    Tetragonal,   // 4/mmm
    Rhombohedral, // -3m
    Hexagonal,    // 6/mmm
    Cubic,        // m-3m
}

pub enum Setting {
    HallNumber(i32),
    /// The setting of the smallest Hall number
    Spglib,
    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral, and origin choice 2 for centrosymmetric space groups
    Standard,
}
