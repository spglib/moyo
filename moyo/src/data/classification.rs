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

/// ===========================================================================
/// Classification based on lattice
/// ===========================================================================
/// Reuse notations for Bravais type of lattices
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum BravaisClass {
    aP,
    mP,
    mC,
    oP,
    oS,
    oF,
    oI,
    tP,
    tI,
    hR,
    hP,
    cP,
    cF,
    cI,
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
