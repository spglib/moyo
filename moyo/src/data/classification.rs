use strum_macros::EnumIter;

/// ===========================================================================
/// Classification based on point group
/// ===========================================================================
/// c.f. Table 3.2.3.2 of ITA(6th)
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum GeometricCrystalClass {
    // Triclinic
    C1, // 1
    Ci, // -1
    // Monoclinic
    C2,  // 2
    C1h, // m
    C2h, // 2/m
    // Orthorhombic
    D2,  // 222
    C2v, // mm2
    D2h, // mmm
    // Tetragonal
    C4,  // 4
    S4,  // -4
    C4h, // 4/m
    D4,  // 422
    C4v, // 4mm
    D2d, // -42m
    D4h, // 4/mmm
    // Trigonal
    C3,  // 3
    C3i, // -3
    D3,  // 32
    C3v, // 3m
    D3d, // -3m
    // Hexagonal
    C6,  // 6
    C3h, // -6
    C6h, // 6/m
    D6,  // 622
    C6v, // 6mm
    D3h, // -6m2
    D6h, // 6/mmm
    // Cubic
    T,  // 23
    Th, // m-3
    O,  // 432
    Td, // -43m
    Oh, // m-3m
}

impl ToString for GeometricCrystalClass {
    fn to_string(&self) -> String {
        match self {
            // Triclinic
            GeometricCrystalClass::C1 => "1".to_string(),
            GeometricCrystalClass::Ci => "-1".to_string(),
            // Monoclinic
            GeometricCrystalClass::C2 => "2".to_string(),
            GeometricCrystalClass::C1h => "m".to_string(),
            GeometricCrystalClass::C2h => "2/m".to_string(),
            // Orthorhombic
            GeometricCrystalClass::D2 => "222".to_string(),
            GeometricCrystalClass::C2v => "mm2".to_string(),
            GeometricCrystalClass::D2h => "mmm".to_string(),
            // Tetragonal
            GeometricCrystalClass::C4 => "4".to_string(),
            GeometricCrystalClass::S4 => "-4".to_string(),
            GeometricCrystalClass::C4h => "4/m".to_string(),
            GeometricCrystalClass::D4 => "422".to_string(),
            GeometricCrystalClass::C4v => "4mm".to_string(),
            GeometricCrystalClass::D2d => "-42m".to_string(),
            GeometricCrystalClass::D4h => "4/mmm".to_string(),
            // Trigonal
            GeometricCrystalClass::C3 => "3".to_string(),
            GeometricCrystalClass::C3i => "-3".to_string(),
            GeometricCrystalClass::D3 => "32".to_string(),
            GeometricCrystalClass::C3v => "3m".to_string(),
            GeometricCrystalClass::D3d => "-3m".to_string(),
            // Hexagonal
            GeometricCrystalClass::C6 => "6".to_string(),
            GeometricCrystalClass::C3h => "-6".to_string(),
            GeometricCrystalClass::C6h => "6/m".to_string(),
            GeometricCrystalClass::D6 => "622".to_string(),
            GeometricCrystalClass::C6v => "6mm".to_string(),
            GeometricCrystalClass::D3h => "-6m2".to_string(),
            GeometricCrystalClass::D6h => "6/mmm".to_string(),
            // Cubic
            GeometricCrystalClass::T => "23".to_string(),
            GeometricCrystalClass::Th => "m-3".to_string(),
            GeometricCrystalClass::O => "432".to_string(),
            GeometricCrystalClass::Td => "-43m".to_string(),
            GeometricCrystalClass::Oh => "m-3m".to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum LaueClass {
    Ci,  // -1
    C2h, // 2/m,
    D2h, // mmm
    C4h, // 4/m
    D4h, // 4/mmm
    C3i, // -3
    D3d, // -3m
    C6h, // 6/m
    D6h, // 6/mmm
    Th,  // m-3
    Oh,  // m-3m
}

impl LaueClass {
    #[allow(dead_code)]
    pub fn from_geometric_crystal_class(geometric_crystal_class: GeometricCrystalClass) -> Self {
        match geometric_crystal_class {
            GeometricCrystalClass::C1 | GeometricCrystalClass::Ci => LaueClass::Ci,
            GeometricCrystalClass::C2 | GeometricCrystalClass::C1h | GeometricCrystalClass::C2h => {
                LaueClass::C2h
            }
            GeometricCrystalClass::D2 | GeometricCrystalClass::C2v | GeometricCrystalClass::D2h => {
                LaueClass::D2h
            }
            GeometricCrystalClass::C4 | GeometricCrystalClass::S4 | GeometricCrystalClass::C4h => {
                LaueClass::C4h
            }
            GeometricCrystalClass::D4
            | GeometricCrystalClass::C4v
            | GeometricCrystalClass::D2d
            | GeometricCrystalClass::D4h => LaueClass::D4h,
            GeometricCrystalClass::C3 | GeometricCrystalClass::C3i => LaueClass::C3i,
            GeometricCrystalClass::D3 | GeometricCrystalClass::C3v | GeometricCrystalClass::D3d => {
                LaueClass::D3d
            }
            GeometricCrystalClass::C6 | GeometricCrystalClass::C3h | GeometricCrystalClass::C6h => {
                LaueClass::C6h
            }
            GeometricCrystalClass::D6
            | GeometricCrystalClass::C6v
            | GeometricCrystalClass::D3h
            | GeometricCrystalClass::D6h => LaueClass::D6h,
            GeometricCrystalClass::T | GeometricCrystalClass::Th => LaueClass::Th,
            GeometricCrystalClass::O | GeometricCrystalClass::Td | GeometricCrystalClass::Oh => {
                LaueClass::Oh
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
}

impl CrystalSystem {
    pub fn from_geometric_crystal_class(geometric_crystal_class: GeometricCrystalClass) -> Self {
        match geometric_crystal_class {
            // Triclinic
            GeometricCrystalClass::C1 | GeometricCrystalClass::Ci => CrystalSystem::Triclinic,
            // Monoclinic
            GeometricCrystalClass::C2 | GeometricCrystalClass::C1h | GeometricCrystalClass::C2h => {
                CrystalSystem::Monoclinic
            }
            // Orthorhombic
            GeometricCrystalClass::D2 | GeometricCrystalClass::C2v | GeometricCrystalClass::D2h => {
                CrystalSystem::Orthorhombic
            }
            // Tetragonal
            GeometricCrystalClass::C4
            | GeometricCrystalClass::S4
            | GeometricCrystalClass::C4h
            | GeometricCrystalClass::D4
            | GeometricCrystalClass::C4v
            | GeometricCrystalClass::D2d
            | GeometricCrystalClass::D4h => CrystalSystem::Tetragonal,
            // Trigonal
            GeometricCrystalClass::C3
            | GeometricCrystalClass::C3i
            | GeometricCrystalClass::D3
            | GeometricCrystalClass::C3v
            | GeometricCrystalClass::D3d => CrystalSystem::Trigonal,
            // Hexagonal
            GeometricCrystalClass::C6
            | GeometricCrystalClass::C3h
            | GeometricCrystalClass::C6h
            | GeometricCrystalClass::D6
            | GeometricCrystalClass::C6v
            | GeometricCrystalClass::D3h
            | GeometricCrystalClass::D6h => CrystalSystem::Hexagonal,
            // Cubic
            GeometricCrystalClass::T
            | GeometricCrystalClass::Th
            | GeometricCrystalClass::O
            | GeometricCrystalClass::Td
            | GeometricCrystalClass::Oh => CrystalSystem::Cubic,
        }
    }
}

impl ToString for CrystalSystem {
    fn to_string(&self) -> String {
        match self {
            CrystalSystem::Triclinic => "Triclinic".to_string(),
            CrystalSystem::Monoclinic => "Monoclinic".to_string(),
            CrystalSystem::Orthorhombic => "Orthorhombic".to_string(),
            CrystalSystem::Tetragonal => "Tetragonal".to_string(),
            CrystalSystem::Trigonal => "Trigonal".to_string(),
            CrystalSystem::Hexagonal => "Hexagonal".to_string(),
            CrystalSystem::Cubic => "Cubic".to_string(),
        }
    }
}

/// ===========================================================================
/// Classification based on lattice
/// ===========================================================================
/// Reuse notations for Bravais type of lattices
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum BravaisClass {
    // Triclinic
    aP,
    // Monoclinic
    mP,
    mC,
    // Orthorhombic
    oP,
    oS,
    oF,
    oI,
    // Tetragonal
    tP,
    tI,
    // Rhombohedral
    hR,
    // Hexagonal
    hP,
    // Cubic
    cP,
    cF,
    cI,
}

impl ToString for BravaisClass {
    fn to_string(&self) -> String {
        match self {
            // Triclinic
            BravaisClass::aP => "aP".to_string(),
            // Monoclinic
            BravaisClass::mP => "mP".to_string(),
            BravaisClass::mC => "mC".to_string(),
            // Orthorhombic
            BravaisClass::oP => "oP".to_string(),
            BravaisClass::oS => "oS".to_string(),
            BravaisClass::oF => "oF".to_string(),
            BravaisClass::oI => "oI".to_string(),
            // Tetragonal
            BravaisClass::tP => "tP".to_string(),
            BravaisClass::tI => "tI".to_string(),
            // Rhombohedral
            BravaisClass::hR => "hR".to_string(),
            // Hexagonal
            BravaisClass::hP => "hP".to_string(),
            // Cubic
            BravaisClass::cP => "cP".to_string(),
            BravaisClass::cF => "cF".to_string(),
            BravaisClass::cI => "cI".to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
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

impl LatticeSystem {
    pub fn from_bravais_class(bravais_class: BravaisClass) -> Self {
        match bravais_class {
            BravaisClass::aP => LatticeSystem::Triclinic,
            BravaisClass::mP | BravaisClass::mC => LatticeSystem::Monoclinic,
            BravaisClass::oP | BravaisClass::oS | BravaisClass::oF | BravaisClass::oI => {
                LatticeSystem::Orthorhombic
            }
            BravaisClass::tP | BravaisClass::tI => LatticeSystem::Tetragonal,
            BravaisClass::hR => LatticeSystem::Rhombohedral,
            BravaisClass::hP => LatticeSystem::Hexagonal,
            BravaisClass::cP | BravaisClass::cF | BravaisClass::cI => LatticeSystem::Cubic,
        }
    }
}

impl ToString for LatticeSystem {
    fn to_string(&self) -> String {
        match self {
            LatticeSystem::Triclinic => "Triclinic".to_string(),
            LatticeSystem::Monoclinic => "Monoclinic".to_string(),
            LatticeSystem::Orthorhombic => "Orthorhombic".to_string(),
            LatticeSystem::Tetragonal => "Tetragonal".to_string(),
            LatticeSystem::Rhombohedral => "Rhombohedral".to_string(),
            LatticeSystem::Hexagonal => "Hexagonal".to_string(),
            LatticeSystem::Cubic => "Cubic".to_string(),
        }
    }
}

/// ===========================================================================
/// Other classification
/// ===========================================================================

#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum CrystalFamily {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Hexagonal,
    Cubic,
}

impl CrystalFamily {
    #[allow(dead_code)]
    pub fn from_crystal_system(crystal_system: CrystalSystem) -> Self {
        match crystal_system {
            CrystalSystem::Triclinic => CrystalFamily::Triclinic,
            CrystalSystem::Monoclinic => CrystalFamily::Monoclinic,
            CrystalSystem::Orthorhombic => CrystalFamily::Orthorhombic,
            CrystalSystem::Tetragonal => CrystalFamily::Tetragonal,
            CrystalSystem::Trigonal | CrystalSystem::Hexagonal => CrystalFamily::Hexagonal,
            CrystalSystem::Cubic => CrystalFamily::Cubic,
        }
    }

    #[allow(dead_code)]
    pub fn from_lattice_system(lattice_system: LatticeSystem) -> Self {
        match lattice_system {
            LatticeSystem::Triclinic => CrystalFamily::Triclinic,
            LatticeSystem::Monoclinic => CrystalFamily::Monoclinic,
            LatticeSystem::Orthorhombic => CrystalFamily::Orthorhombic,
            LatticeSystem::Tetragonal => CrystalFamily::Tetragonal,
            LatticeSystem::Rhombohedral | LatticeSystem::Hexagonal => CrystalFamily::Hexagonal,
            LatticeSystem::Cubic => CrystalFamily::Cubic,
        }
    }
}

impl ToString for CrystalFamily {
    fn to_string(&self) -> String {
        match self {
            CrystalFamily::Triclinic => "Triclinic".to_string(),
            CrystalFamily::Monoclinic => "Monoclinic".to_string(),
            CrystalFamily::Orthorhombic => "Orthorhombic".to_string(),
            CrystalFamily::Tetragonal => "Tetragonal".to_string(),
            CrystalFamily::Hexagonal => "Hexagonal".to_string(),
            CrystalFamily::Cubic => "Cubic".to_string(),
        }
    }
}
