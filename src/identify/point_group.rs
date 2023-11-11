use crate::base::lattice::Lattice;
use crate::base::operation::Rotation;

/// Crystallographic point group
#[derive(Debug)]
pub struct PointGroup {
    pub lattice: Lattice,
    pub rotations: Vec<Rotation>,
}

impl PointGroup {
    pub fn new(lattice: Lattice, rotations: Vec<Rotation>) -> Self {
        Self { lattice, rotations }
    }
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

enum RotationType {
    Rotation1,      // 1
    Rotation2,      // 2
    Rotation3,      // 3
    Rotation4,      // 4
    Rotation6,      // 6
    RotoInversion1, // -1 = S2
    RotoInversion2, // -2 = m = S1
    RotoInversion3, // -3 = S6^-1
    RotoInversion4, // -4 = S4^-1
    RotoInversion6, // -6 = S3^-1
}

fn identify_rotation_type(rotation: &Rotation) -> Option<RotationType> {
    let tr = rotation.trace();
    let det = rotation.map(|e| e as f64).determinant().round() as i32;

    match (tr, det) {
        (3, 1) => Some(RotationType::Rotation1),
        (-1, 1) => Some(RotationType::Rotation2),
        (0, 1) => Some(RotationType::Rotation3),
        (1, 1) => Some(RotationType::Rotation4),
        (2, 1) => Some(RotationType::Rotation6),
        (-3, -1) => Some(RotationType::RotoInversion1),
        (1, -1) => Some(RotationType::RotoInversion2),
        (0, -1) => Some(RotationType::RotoInversion3),
        (-1, -1) => Some(RotationType::RotoInversion4),
        (-2, -1) => Some(RotationType::RotoInversion6),
        _ => None,
    }
}

/// Use look up table in Table 6 of https://arxiv.org/pdf/1808.01590.pdf
pub fn identify_geometric_crystal_class(point_group: &PointGroup) -> Option<GeometricCrystalClass> {
    // count RotationTypes in point_group
    let mut rotation_types = [0; 10];
    for rotation in &point_group.rotations {
        if let Some(rotation_type) = identify_rotation_type(rotation) {
            match rotation_type {
                RotationType::RotoInversion6 => rotation_types[0] += 1,
                RotationType::RotoInversion4 => rotation_types[1] += 1,
                RotationType::RotoInversion3 => rotation_types[2] += 1,
                RotationType::RotoInversion2 => rotation_types[3] += 1,
                RotationType::RotoInversion1 => rotation_types[4] += 1,
                RotationType::Rotation1 => rotation_types[5] += 1,
                RotationType::Rotation2 => rotation_types[6] += 1,
                RotationType::Rotation3 => rotation_types[7] += 1,
                RotationType::Rotation4 => rotation_types[8] += 1,
                RotationType::Rotation6 => rotation_types[9] += 1,
            }
        } else {
            return None;
        }
    }
    match rotation_types {
        // Triclinic
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::C1),
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::Ci),
        // Monoclinic
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2),
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0] => Some(GeometricCrystalClass::C1h),
        [0, 0, 0, 1, 1, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2h),
        // Orthorhombic
        [0, 0, 0, 0, 0, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2),
        [0, 0, 0, 2, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::C2v),
        [0, 0, 0, 3, 1, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2h),
        // Tetragonal
        [0, 0, 0, 0, 0, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4),
        [0, 2, 0, 0, 0, 1, 1, 0, 0, 0] => Some(GeometricCrystalClass::S4),
        [0, 2, 0, 1, 1, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4h),
        [0, 0, 0, 0, 0, 1, 5, 0, 2, 0] => Some(GeometricCrystalClass::D4),
        [0, 0, 0, 4, 0, 1, 1, 0, 2, 0] => Some(GeometricCrystalClass::C4v),
        [0, 2, 0, 2, 0, 1, 3, 0, 0, 0] => Some(GeometricCrystalClass::D2d),
        [0, 2, 0, 5, 1, 1, 5, 0, 2, 0] => Some(GeometricCrystalClass::D4h),
        // Trigonal
        [0, 0, 0, 0, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3),
        [0, 0, 2, 0, 1, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3i),
        [0, 0, 0, 0, 0, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3),
        [0, 0, 0, 3, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3v),
        [0, 0, 2, 3, 1, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3d),
        // Hexagonal
        [0, 0, 0, 0, 0, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6),
        [2, 0, 0, 1, 0, 1, 0, 2, 0, 0] => Some(GeometricCrystalClass::C3h),
        [2, 0, 2, 1, 1, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6h),
        [0, 0, 0, 0, 0, 1, 7, 2, 0, 2] => Some(GeometricCrystalClass::D6),
        [0, 0, 0, 6, 0, 1, 1, 2, 0, 2] => Some(GeometricCrystalClass::C6v),
        [2, 0, 0, 4, 0, 1, 3, 2, 0, 0] => Some(GeometricCrystalClass::D3h),
        [2, 0, 2, 7, 1, 1, 7, 2, 0, 2] => Some(GeometricCrystalClass::D6h),
        // Cubic
        [0, 0, 0, 0, 0, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::T),
        [0, 0, 8, 3, 1, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::Th),
        [0, 0, 0, 0, 0, 1, 9, 8, 6, 0] => Some(GeometricCrystalClass::O),
        [0, 6, 0, 6, 0, 1, 3, 8, 0, 0] => Some(GeometricCrystalClass::Td),
        [0, 6, 8, 9, 1, 1, 9, 8, 6, 0] => Some(GeometricCrystalClass::Oh),
        // Unknown
        _ => None,
    }
}

#[cfg(test)]
mod tests {}
