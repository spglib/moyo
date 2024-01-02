use rstest::rstest;

use moyo::base::cell::Cell;
use moyo::base::tolerance::AngleTolerance;
use moyo::search::primitive_cell::PrimitiveCellSearch;
use moyo::search::symmetry_search::SymmetrySearch;

mod fixtures;
use fixtures::rutile;

#[rstest]
fn test_conventional_cell_rutile(rutile: Cell) {
    let cell = rutile;

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let primitive_cell_result = PrimitiveCellSearch::new(&cell, symprec).unwrap();
    let primitive_cell = primitive_cell_result.primitive_cell;

    let _ = SymmetrySearch::new(&primitive_cell, symprec, angle_tolerance).unwrap();
}
