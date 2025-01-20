use super::primitive_cell::{PrimitiveCell, PrimitiveMagneticCell};
use super::primitive_symmetry_search::{PrimitiveMagneticSymmetrySearch, PrimitiveSymmetrySearch};
use crate::base::{
    AngleTolerance, Cell, MagneticCell, MagneticMoment, MagneticSymmetryTolerances, MoyoError,
    RotationMagneticMomentAction, SymmetryTolerances, ToleranceHandler,
};

use log::debug;

const MAX_SYMMETRY_SEARCH_TRIALS: usize = 16;
const MAX_TOLERANCE_HANDLER_TRIALS: usize = 4;

pub fn iterative_symmetry_search(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<(PrimitiveCell, PrimitiveSymmetrySearch, f64, AngleTolerance), MoyoError> {
    let mut tolerances = SymmetryTolerances {
        symprec,
        angle_tolerance,
    };

    for _ in 0..MAX_TOLERANCE_HANDLER_TRIALS {
        let mut tolerance_handler = ToleranceHandler::new(tolerances);

        for _ in 0..MAX_SYMMETRY_SEARCH_TRIALS {
            match PrimitiveCell::new(cell, tolerance_handler.tolerances.symprec) {
                Ok(prim_cell) => {
                    match PrimitiveSymmetrySearch::new(
                        &prim_cell.cell,
                        tolerance_handler.tolerances.symprec,
                        tolerance_handler.tolerances.angle_tolerance,
                    ) {
                        Ok(symmetry_search) => {
                            return Ok((
                                prim_cell,
                                symmetry_search,
                                tolerance_handler.tolerances.symprec,
                                tolerance_handler.tolerances.angle_tolerance,
                            ));
                        }
                        Err(err) => tolerance_handler.update(err),
                    }
                }
                Err(err) => tolerance_handler.update(err),
            }
        }

        // When the maximum number of symmetry search trials is reached, restart ToleranceHandler to try larger strides
        tolerances = tolerance_handler.tolerances.clone();
        debug!("Restart ToleranceHandler with {:?}", tolerances);
    }
    debug!("Reach the maximum number of symmetry search trials");
    Err(MoyoError::PrimitiveSymmetrySearchError)
}

pub fn iterative_magnetic_symmetry_search<M: MagneticMoment>(
    magnetic_cell: &MagneticCell<M>,
    symprec: f64,
    angle_tolerance: AngleTolerance,
    mag_symprec: Option<f64>,
    action: RotationMagneticMomentAction,
) -> Result<
    (
        PrimitiveMagneticCell<M>,
        PrimitiveMagneticSymmetrySearch,
        f64,
        AngleTolerance,
        f64,
    ),
    MoyoError,
> {
    let mut tolerances = MagneticSymmetryTolerances {
        symprec,
        angle_tolerance,
        mag_symprec: mag_symprec.unwrap_or(symprec),
    };

    for _ in 0..MAX_TOLERANCE_HANDLER_TRIALS {
        let mut tolerance_handler = ToleranceHandler::new(tolerances);

        for _ in 0..MAX_SYMMETRY_SEARCH_TRIALS {
            match PrimitiveMagneticCell::new(
                magnetic_cell,
                tolerance_handler.tolerances.symprec,
                tolerance_handler.tolerances.mag_symprec,
            ) {
                Ok(prim_mag_cell) => {
                    match PrimitiveMagneticSymmetrySearch::new(
                        &prim_mag_cell.magnetic_cell,
                        tolerance_handler.tolerances.symprec,
                        tolerance_handler.tolerances.angle_tolerance,
                        tolerance_handler.tolerances.mag_symprec,
                        action,
                    ) {
                        Ok(magnetic_symmetry_search) => {
                            return Ok((
                                prim_mag_cell,
                                magnetic_symmetry_search,
                                tolerance_handler.tolerances.symprec,
                                tolerance_handler.tolerances.angle_tolerance,
                                tolerance_handler.tolerances.mag_symprec,
                            ));
                        }
                        Err(err) => tolerance_handler.update(err),
                    }
                }
                Err(err) => tolerance_handler.update(err),
            }
        }

        // When the maximum number of symmetry search trials is reached, restart ToleranceHandler to try larger strides
        tolerances = tolerance_handler.tolerances.clone();
        debug!("Restart ToleranceHandler with {:?}", tolerances);
    }
    debug!("Reach the maximum number of symmetry search trials");
    Err(MoyoError::PrimitiveMagneticSymmetrySearchError)
}
