use crate::base::operation::{AbstractOperations, Operations};

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Setting {
    HallNumber(i32),
    /// The setting of the smallest Hall number
    Spglib,
    /// Unique axis b, cell choice 1 for monoclinic, hexagonal axes for rhombohedral, and origin choice 2 for centrosymmetric space groups
    Standard,
}

pub struct SpaceGroupRepresentative {
    pub generators: AbstractOperations,
}

pub struct SpaceGroup {
    pub operations: Operations,
}
