use rstest::*;

use nalgebra::{matrix, Vector3};

use moyo::base::cell::Cell;
use moyo::base::lattice::Lattice;

#[fixture]
pub fn primitive_fcc() -> Cell {
    let lattice = Lattice::new(matrix![
            0.0, 0.5, 0.5;
            0.5, 0.0, 0.5;
            0.5, 0.5, 0.0;
    ]);
    let positions = vec![Vector3::new(0.0, 0.0, 0.0)];
    let numbers = vec![0];

    Cell::new(lattice, positions, numbers)
}

#[fixture]
/// Rutile, P4_2/mnm (No. 136)
pub fn rutile() -> Cell {
    let a = 4.603;
    let c = 2.969;
    let lattice = Lattice::new(matrix![
        a, 0.0, 0.0;
        0.0, a, 0.0;
        0.0, 0.0, c;
    ]);
    let x_4f = 0.3046;
    let positions = vec![
        Vector3::new(0.0, 0.0, 0.0),                // Ti(2a)
        Vector3::new(0.5, 0.5, 0.5),                // Ti(2a)
        Vector3::new(x_4f, x_4f, 0.0),              // O(4f)
        Vector3::new(-x_4f, -x_4f, 0.0),            // O(4f)
        Vector3::new(-x_4f + 0.5, x_4f + 0.5, 0.5), // O(4f)
        Vector3::new(x_4f + 0.5, -x_4f + 0.5, 0.5), // O(4f)
    ];
    let numbers = vec![0, 0, 1, 1, 1, 1];

    Cell::new(lattice, positions, numbers)
}
