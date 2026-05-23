use super::layer::{LayerPrimitiveCell, LayerPrimitiveSymmetrySearch};
use super::primitive_cell::{PrimitiveCell, PrimitiveMagneticCell};
use super::primitive_symmetry_search::{PrimitiveMagneticSymmetrySearch, PrimitiveSymmetrySearch};
use crate::base::{
    AngleTolerance, Cell, LayerCell, MagneticCell, MagneticMoment, MagneticSymmetryTolerances,
    MoyoError, RotationMagneticMomentAction, SymmetryTolerances, ToleranceHandler,
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

/// Layer-group counterpart of [`iterative_symmetry_search`].
///
/// `LayerCell::new` validates the structural layer-cell contract (`c`
/// perpendicular to `a, b`); the errors it raises (e.g.
/// [`MoyoError::AperiodicAxisNotOrthogonal`]) describe the input geometry, not
/// the tolerance, so it sits outside the retry loop. Only
/// [`LayerPrimitiveCell::new`] and [`LayerPrimitiveSymmetrySearch::new`] are
/// retried under tightening tolerances, matching the bulk path's choice to
/// wrap only the tolerance-sensitive symmetry-search stages.
pub fn iterative_layer_symmetry_search(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
) -> Result<
    (
        LayerPrimitiveCell,
        LayerPrimitiveSymmetrySearch,
        f64,
        AngleTolerance,
    ),
    MoyoError,
> {
    let layer_cell = LayerCell::new(cell.clone(), symprec, angle_tolerance)?;

    let mut tolerances = SymmetryTolerances {
        symprec,
        angle_tolerance,
    };

    for _ in 0..MAX_TOLERANCE_HANDLER_TRIALS {
        let mut tolerance_handler = ToleranceHandler::new(tolerances);

        for _ in 0..MAX_SYMMETRY_SEARCH_TRIALS {
            match LayerPrimitiveCell::new(&layer_cell, tolerance_handler.tolerances.symprec) {
                Ok(prim_layer) => {
                    match LayerPrimitiveSymmetrySearch::new(
                        &prim_layer.layer_cell,
                        tolerance_handler.tolerances.symprec,
                        tolerance_handler.tolerances.angle_tolerance,
                    ) {
                        Ok(symmetry_search) => {
                            return Ok((
                                prim_layer,
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

#[cfg(test)]
mod tests {
    use nalgebra::{Vector3, matrix};

    use super::iterative_layer_symmetry_search;
    use crate::base::{AngleTolerance, Cell, Lattice, LayerCell, MoyoError};
    use crate::search::layer::LayerPrimitiveCell;

    fn unit_square_p1_cell() -> Cell {
        // p1 square layer, unit in-plane basis. `LayerPrimitiveCell::new`
        // raises `TooLargeToleranceError` whenever `2 * symprec > 0.5`.
        Cell::new(
            Lattice::new(matrix![
                1.0, 0.0, 0.0;
                0.0, 1.0, 0.0;
                0.0, 0.0, 5.0;
            ]),
            vec![Vector3::new(0.3, 0.4, 0.2)],
            vec![1],
        )
    }

    #[test]
    fn test_iterative_layer_symmetry_search_passthrough_when_symprec_is_already_safe() {
        // With a tight symprec the inner search succeeds on the first try, so
        // the returned tolerances must equal the inputs.
        let cell = unit_square_p1_cell();
        let symprec = 1e-4;
        let angle = AngleTolerance::Default;

        let (_prim, _search, out_symprec, out_angle) =
            iterative_layer_symmetry_search(&cell, symprec, angle).unwrap();
        assert_eq!(out_symprec, symprec);
        match (out_angle, angle) {
            (AngleTolerance::Default, AngleTolerance::Default) => {}
            other => panic!("angle tolerance changed unexpectedly: {:?}", other),
        }
    }

    #[test]
    fn test_iterative_layer_symmetry_search_retries_after_too_large_tolerance() {
        // For this cell, `LayerPrimitiveCell::new` errors directly at symprec
        // = 0.3 (since 2 * 0.3 > 0.5). Confirm that error in the direct call
        // first, then confirm the iterative wrapper recovers by reducing
        // symprec and returning a tightened value.
        let cell = unit_square_p1_cell();
        let layer = LayerCell::new(cell.clone(), 0.3, AngleTolerance::Default).unwrap();
        assert!(matches!(
            LayerPrimitiveCell::new(&layer, 0.3),
            Err(MoyoError::TooLargeToleranceError),
        ));

        let (_prim, _search, out_symprec, _out_angle) =
            iterative_layer_symmetry_search(&cell, 0.3, AngleTolerance::Default).unwrap();
        assert!(
            out_symprec < 0.3,
            "expected the retry to tighten symprec below the input; got {}",
            out_symprec,
        );
    }
}
