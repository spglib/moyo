use thiserror::Error;

#[derive(Error, Debug)]
pub enum MoyoError {
    #[error("Minkowski reduction failed")]
    MinkowskiReductionError,
    #[error("Too small symprec")]
    TooSmallSymprecError,
    #[error("Primitive cell search failed")]
    PrimitiveCellSearchError,
    #[error("Bravais group search failed")]
    BravaisGroupSearchError,
}
