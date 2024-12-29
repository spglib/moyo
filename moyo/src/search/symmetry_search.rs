use super::primitive_cell::PrimitiveCell;
use super::primitive_symmetry_search::PrimitiveSymmetrySearch;
use crate::base::{AngleTolerance, Cell, MoyoError, ToleranceHandler};

use log::debug;

const MAX_SYMMETRY_SEARCH_TRIALS: usize = 16;
const MAX_TOLERANCE_HANDLER_TRIALS: usize = 4;

pub fn iterative_symmetry_search(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<(PrimitiveCell, PrimitiveSymmetrySearch, f64, AngleTolerance), MoyoError> {
    let mut current_symprec = symprec;
    let mut current_angle_tolerance = angle_tolerance;

    for _ in 0..MAX_TOLERANCE_HANDLER_TRIALS {
        let mut tolerance_handler = ToleranceHandler::new(current_symprec, current_angle_tolerance);

        for _ in 0..MAX_SYMMETRY_SEARCH_TRIALS {
            match PrimitiveCell::new(cell, tolerance_handler.symprec) {
                Ok(prim_cell) => {
                    match PrimitiveSymmetrySearch::new(
                        &prim_cell.cell,
                        tolerance_handler.symprec,
                        tolerance_handler.angle_tolerance,
                    ) {
                        Ok(symmetry_search) => {
                            return Ok((
                                prim_cell,
                                symmetry_search,
                                tolerance_handler.symprec,
                                tolerance_handler.angle_tolerance,
                            ));
                        }
                        Err(err) => tolerance_handler.update(err),
                    }
                }
                Err(err) => tolerance_handler.update(err),
            }
        }

        // When the maximum number of symmetry search trials is reached, restart ToleranceHandler to try larger strides
        current_symprec = tolerance_handler.symprec;
        current_angle_tolerance = tolerance_handler.angle_tolerance;
        debug!(
            "Restart ToleranceHandler with symprec={}, angle_tolerance={:?}",
            current_symprec, current_angle_tolerance
        );
    }
    debug!("Reach the maximum number of symmetry search trials");
    Err(MoyoError::PrimitiveSymmetrySearchError)
}
