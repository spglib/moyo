use thiserror::Error;

#[derive(Error, Debug, PartialEq, Clone, Copy)]
/// Error types for the **moyo** library
pub enum MoyoError {
    // Lattice reduction errors
    #[error("Minkowski reduction failed")]
    MinkowskiReductionError,
    #[error("Niggli reduction failed")]
    NiggliReductionError,
    #[error("Delaunay reduction failed")]
    DelaunayReductionError,
    // Tolerance related errors
    #[error("Too small tolerance")]
    TooSmallToleranceError,
    #[error("Too large tolerance")]
    TooLargeToleranceError,
    // Symmetry search errors
    #[error("Primitive cell search failed")]
    PrimitiveCellError,
    #[error("Primitive symmetry search failed")]
    PrimitiveSymmetrySearchError,
    #[error("Primitive magnetic symmetry search failed")]
    PrimitiveMagneticSymmetrySearchError,
    #[error("Bravais group search failed")]
    BravaisGroupSearchError,
    // Identification errors
    #[error(
        "Geometric crystal class identification failed. Given rotations might not form a group or contain duplicates."
    )]
    GeometricCrystalClassIdentificationError,
    #[error("Arithmetic crystal class identification failed")]
    ArithmeticCrystalClassIdentificationError,
    #[error("Space group type identification failed")]
    SpaceGroupTypeIdentificationError,
    #[error("Layer group type identification failed")]
    LayerGroupTypeIdentificationError,
    #[error("Construct type identification failed")]
    ConstructTypeIdentificationError,
    #[error("Magnetic space group type identification failed")]
    MagneticSpaceGroupTypeIdentificationError,
    // Standardization errors
    #[error("Standardization failed")]
    StandardizationError,
    #[error("Magnetic standardization failed")]
    MagneticStandardizationError,
    #[error("Wyckoff position assignment failed")]
    WyckoffPositionAssignmentError,
    // Dataset errors
    #[error("Hall symbol parsing failed")]
    HallSymbolParsingError,
    #[error("Unknown hall_number")]
    UnknownHallNumberError,
    #[error("Unknown number")]
    UnknownNumberError,
    #[error("Unknown arithmetic_number")]
    UnknownArithmeticNumberError,
    #[error("Unknown uni_number")]
    UnknownUNINumberError,
    // Layer-group input validation errors
    #[error("Lattice basis vector has zero norm")]
    DegenerateLattice,
    #[error(
        "Aperiodic axis c is not perpendicular to the in-plane axes within tolerance: dev(c,a)={dev_ca:.6} rad, dev(c,b)={dev_cb:.6} rad"
    )]
    AperiodicAxisNotOrthogonal { dev_ca: f64, dev_cb: f64 },
    #[error(
        "In-plane axes a, b are not in the xy-plane within tolerance: dev(a,xy)={dev_az:.6} rad, dev(b,xy)={dev_bz:.6} rad"
    )]
    InPlaneAxesNotInXY { dev_az: f64, dev_bz: f64 },
}
