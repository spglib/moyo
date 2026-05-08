use strum_macros::EnumIter;

use super::layer_centering::LayerCentering;

/// 2D lattice systems applicable to layer groups.
/// Reference: Fu et al. 2024, Table 3.
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum LayerLatticeSystem {
    Oblique,     // 2/m
    Rectangular, // mmm
    Square,      // 4/mmm
    Hexagonal,   // 6/mmm
}

#[allow(clippy::to_string_trait_impl)]
impl ToString for LayerLatticeSystem {
    fn to_string(&self) -> String {
        match self {
            LayerLatticeSystem::Oblique => "Oblique".to_string(),
            LayerLatticeSystem::Rectangular => "Rectangular".to_string(),
            LayerLatticeSystem::Square => "Square".to_string(),
            LayerLatticeSystem::Hexagonal => "Hexagonal".to_string(),
        }
    }
}

/// Bravais types of 2D lattices for layer groups.
/// Reference: Fu et al. 2024, Table 3.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum LayerBravaisClass {
    mp, // Oblique primitive
    op, // Rectangular primitive
    oc, // Rectangular centered
    tp, // Square primitive
    hp, // Hexagonal primitive
}

#[allow(clippy::to_string_trait_impl)]
impl ToString for LayerBravaisClass {
    fn to_string(&self) -> String {
        match self {
            LayerBravaisClass::mp => "mp".to_string(),
            LayerBravaisClass::op => "op".to_string(),
            LayerBravaisClass::oc => "oc".to_string(),
            LayerBravaisClass::tp => "tp".to_string(),
            LayerBravaisClass::hp => "hp".to_string(),
        }
    }
}

impl LayerLatticeSystem {
    pub fn from_layer_bravais_class(layer_bravais_class: LayerBravaisClass) -> Self {
        match layer_bravais_class {
            LayerBravaisClass::mp => LayerLatticeSystem::Oblique,
            LayerBravaisClass::op | LayerBravaisClass::oc => LayerLatticeSystem::Rectangular,
            LayerBravaisClass::tp => LayerLatticeSystem::Square,
            LayerBravaisClass::hp => LayerLatticeSystem::Hexagonal,
        }
    }
}

impl LayerCentering {
    pub fn from_layer_bravais_class(layer_bravais_class: LayerBravaisClass) -> Self {
        match layer_bravais_class {
            LayerBravaisClass::mp
            | LayerBravaisClass::op
            | LayerBravaisClass::tp
            | LayerBravaisClass::hp => LayerCentering::P,
            LayerBravaisClass::oc => LayerCentering::C,
        }
    }
}

/// Crystal systems applicable to layer groups.
/// Reference: Fu et al. 2024, Table 2 (point-group classification of layer groups).
/// The oblique-vs-rectangular distinction for monoclinic is a *lattice-system*
/// distinction (Table 3), not a crystal-system one.
#[derive(Debug, Clone, Copy, PartialEq, EnumIter)]
pub enum LayerCrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
}

#[allow(clippy::to_string_trait_impl)]
impl ToString for LayerCrystalSystem {
    fn to_string(&self) -> String {
        match self {
            LayerCrystalSystem::Triclinic => "Triclinic".to_string(),
            LayerCrystalSystem::Monoclinic => "Monoclinic".to_string(),
            LayerCrystalSystem::Orthorhombic => "Orthorhombic".to_string(),
            LayerCrystalSystem::Tetragonal => "Tetragonal".to_string(),
            LayerCrystalSystem::Trigonal => "Trigonal".to_string(),
            LayerCrystalSystem::Hexagonal => "Hexagonal".to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use strum::IntoEnumIterator;

    use super::*;

    #[test]
    fn test_lattice_system_from_bravais() {
        assert_eq!(
            LayerLatticeSystem::from_layer_bravais_class(LayerBravaisClass::mp),
            LayerLatticeSystem::Oblique
        );
        assert_eq!(
            LayerLatticeSystem::from_layer_bravais_class(LayerBravaisClass::op),
            LayerLatticeSystem::Rectangular
        );
        assert_eq!(
            LayerLatticeSystem::from_layer_bravais_class(LayerBravaisClass::oc),
            LayerLatticeSystem::Rectangular
        );
        assert_eq!(
            LayerLatticeSystem::from_layer_bravais_class(LayerBravaisClass::tp),
            LayerLatticeSystem::Square
        );
        assert_eq!(
            LayerLatticeSystem::from_layer_bravais_class(LayerBravaisClass::hp),
            LayerLatticeSystem::Hexagonal
        );
    }

    #[test]
    fn test_centering_from_bravais() {
        assert_eq!(
            LayerCentering::from_layer_bravais_class(LayerBravaisClass::oc),
            LayerCentering::C
        );
        for bc in [
            LayerBravaisClass::mp,
            LayerBravaisClass::op,
            LayerBravaisClass::tp,
            LayerBravaisClass::hp,
        ] {
            assert_eq!(
                LayerCentering::from_layer_bravais_class(bc),
                LayerCentering::P
            );
        }
    }

    #[test]
    fn test_to_string_round_trip() {
        for ls in LayerLatticeSystem::iter() {
            assert!(!ls.to_string().is_empty());
        }
        for bc in LayerBravaisClass::iter() {
            assert!(!bc.to_string().is_empty());
        }
        for cs in LayerCrystalSystem::iter() {
            assert!(!cs.to_string().is_empty());
        }
    }
}
