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
    #[error("Symmetry search failed")]
    SymmetrySearchError,
    #[error("Arithmetic crystal class identification failed")]
    ArithmeticCrystalClassIdentificationError,
    #[error("Space group type identification failed")]
    SpaceGroupTypeIdentificationError,
}
