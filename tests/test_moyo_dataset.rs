use nalgebra::{matrix, vector, Matrix3, Vector3};

use moyo::base::cell::Cell;
use moyo::base::lattice::Lattice;
use moyo::base::tolerance::AngleTolerance;
use moyo::data::setting::Setting;
use moyo::MoyoDataset;

#[test]
fn test_with_fcc() {
    let lattice = Lattice::new(Matrix3::identity());
    let positions = vec![
        vector![0.0, 0.0, 0.0],
        vector![0.0, 0.5, 0.5],
        vector![0.5, 0.0, 0.5],
        vector![0.5, 0.5, 0.0],
    ];
    let numbers = vec![0, 0, 0, 0];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let dataset =
        MoyoDataset::new(&cell, symprec, AngleTolerance::Default, Setting::Spglib).unwrap();

    assert_eq!(dataset.number, 225); // Fm-3m
    assert_eq!(dataset.hall_number, 523);
    assert_eq!(dataset.num_operations(), 48 * 4);
    assert_eq!(dataset.orbits, vec![0, 0, 0, 0]);
}

#[test]
fn test_with_rutile() {
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
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let dataset =
        MoyoDataset::new(&cell, symprec, AngleTolerance::Default, Setting::Spglib).unwrap();

    assert_eq!(dataset.number, 136); // P4_2/mnm
    assert_eq!(dataset.hall_number, 419);
    assert_eq!(dataset.num_operations(), 16);
}
