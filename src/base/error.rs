use thiserror::Error;

#[derive(Error, Debug)]
pub enum MoyoError {
    #[error("Minkowski reduction failed")]
    MinkowskiReductionError,
    #[error("Niggli reduction failed")]
    NiggliReductionError,
    #[error("Delaunay reduction failed")]
    DelaunayReductionError,
    #[error("Too small symprec")]
    TooSmallSymprecError,
    #[error("Primitive cell search failed")]
    PrimitiveCellError,
    #[error("Bravais group search failed")]
    BravaisGroupSearchError,
    #[error("Primitive symmetry search failed")]
    PrimitiveSymmetrySearchError,
    #[error("Arithmetic crystal class identification failed")]
    ArithmeticCrystalClassIdentificationError,
    #[error("Space group type identification failed")]
    SpaceGroupTypeIdentificationError,
    #[error("Standardization failed")]
    StandardizationError,
    #[error("Wyckoff position assignment failed")]
    WyckoffPositionAssignmentError,
}
