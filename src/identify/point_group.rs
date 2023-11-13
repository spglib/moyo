use crate::base::lattice::Lattice;
use crate::base::operation::{Rotation, RotationType};
use crate::data::classification::GeometricCrystalClass;

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
