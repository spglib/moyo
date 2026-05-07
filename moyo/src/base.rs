mod action;
mod cell;
mod error;
mod lattice;
mod layer_cell;
mod layer_lattice;
mod magnetic_cell;
mod operation;
mod permutation;
mod tolerance;
mod transformation;

pub use action::RotationMagneticMomentAction;
pub use cell::{AtomicSpecie, Cell, Position, validate_aperiodic_axis};
pub use error::MoyoError;
pub use lattice::Lattice;
pub use layer_cell::LayerCell;
pub use layer_lattice::LayerLattice;
pub use magnetic_cell::{Collinear, MagneticCell, MagneticMoment, NonCollinear};
pub use operation::{
    MagneticOperation, MagneticOperations, Operation, Operations, Rotation, Rotations,
    TimeReversal, Translation,
};
pub use permutation::Permutation;
pub use tolerance::{AngleTolerance, is_angle_within_tolerance};
pub use transformation::{Linear, OriginShift, Transformation, UnimodularTransformation};

pub(super) use cell::orbits_from_permutations;
pub(super) use operation::project_rotations;
#[allow(unused_imports)]
pub(super) use operation::traverse;
pub(super) use tolerance::{EPS, MagneticSymmetryTolerances, SymmetryTolerances, ToleranceHandler};
pub(super) use transformation::UnimodularLinear;
