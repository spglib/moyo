use rstest::rstest;

use moyo::base::cell::Cell;
use moyo::base::tolerance::AngleTolerance;
use moyo::search::primitive_cell::search_primitive_cell;
use moyo::search::symmetry_search::search_symmetry_operations_from_primitive;

mod fixtures;
use fixtures::{primitive_fcc, rutile};

#[rstest]
fn test_search_symmetry_operations_fcc(primitive_fcc: Cell) {
    let primitive_cell = primitive_fcc;

    let symprec = 1e-4;
    let result = search_symmetry_operations_from_primitive(
        &primitive_cell,
        symprec,
        AngleTolerance::Default,
    )
    .unwrap();
    assert_eq!(result.operations.num_operations(), 48);
}

#[rstest]
fn test_search_symmetry_operations_rutile(rutile: Cell) {
    let cell = rutile;

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let primitive_cell_result = search_primitive_cell(&cell, symprec).unwrap();
    let primitive_cell = primitive_cell_result.primitive_cell;

    let result =
        search_symmetry_operations_from_primitive(&primitive_cell, symprec, angle_tolerance)
            .unwrap();
    assert_eq!(result.bravais_group.len(), 16);
    assert_eq!(result.operations.num_operations(), 16);
}
