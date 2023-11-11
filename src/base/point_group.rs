use super::lattice::Lattice;
use super::operation::Rotation;

/// Crystallographic point group
#[derive(Debug)]
pub struct PointGroup {
    pub lattice: Lattice,
    pub rotations: Vec<Rotation>,
}

pub enum LaueClass {
    Ci,  // -1
    C2h, // 2/m,
    D2h, // mmm
    C4h, // 4/m
    D4h, // 4/mmm
    C6h, // 6/m
    D6h, // 6/mmm
    Th,  // m-3
    Oh,  // m-3m
}

pub enum CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
}

/// c.f. Table 3.2.3.2 of ITA(6th)
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

impl CrystalSystem {
    pub fn from_geometric_crystal_class(geometric_crystal_class: &GeometricCrystalClass) -> Self {
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

impl PointGroup {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>) -> Self {
        Self { lattice, rotations }
    }
}
