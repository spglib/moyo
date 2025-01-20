use thiserror::Error;

#[derive(Error, Debug, PartialEq, Eq, Clone, Copy)]
/// Error types for the **moyo** library
pub enum MoyoError {
    #[error("Minkowski reduction failed")]
    MinkowskiReductionError,
    #[error("Niggli reduction failed")]
    NiggliReductionError,
    #[error("Delaunay reduction failed")]
    DelaunayReductionError,
    #[error("Too small tolerance")]
    TooSmallToleranceError,
    #[error("Too large tolerance")]
    TooLargeToleranceError,
    #[error("Primitive cell search failed")]
    PrimitiveCellError,
    #[error("Primitive symmetry search failed")]
    PrimitiveSymmetrySearchError,
    #[error("Primitive magnetic symmetry search failed")]
    PrimitiveMagneticSymmetrySearchError,
    #[error("Bravais group search failed")]
    BravaisGroupSearchError,
    #[error("Geometric crystal class identification failed")]
    GeometricCrystalClassIdentificationError,
    #[error("Arithmetic crystal class identification failed")]
    ArithmeticCrystalClassIdentificationError,
    #[error("Space group type identification failed")]
    SpaceGroupTypeIdentificationError,
    #[error("Construct type identification failed")]
    ConstructTypeIdentificationError,
    #[error("Magnetic space group type identification failed")]
    MagneticSpaceGroupTypeIdentificationError,
    #[error("Standardization failed")]
    StandardizationError,
    #[error("Magnetic standardization failed")]
    MagneticStandardizationError,
    #[error("Wyckoff position assignment failed")]
    WyckoffPositionAssignmentError,
    #[error("Hall symbol parsing failed")]
    HallSymbolParsingError,
    #[error("Unknown hall_number")]
    UnknownHallNumberError,
    #[error("Unknown number")]
    UnknownNumberError,
}
