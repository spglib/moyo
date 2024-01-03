use nalgebra::{matrix, vector, Matrix3, Vector3};

use moyo::base::cell::Cell;
use moyo::base::lattice::Lattice;
use moyo::base::operation::AbstractOperations;
use moyo::base::tolerance::AngleTolerance;
use moyo::data::setting::Setting;
use moyo::MoyoDataset;

/// O(num_atoms^2)
fn check_operations(cell: &Cell, operations: &AbstractOperations) -> bool {
    // Check uniqueness
    let num_operations = operations.num_operations();
    for i in 0..num_operations {
        for j in i + 1..num_operations {
            if operations.rotations[i] != operations.rotations[j] {
                continue;
            }
            let mut diff = operations.translations[i] - operations.translations[j];
            diff -= diff.map(|x| x.round());
            if diff.iter().all(|x| x.abs() < 1e-4) {
                return false;
            }
        }
    }

    // Check if overlap
    for (rotation, translation) in operations
        .rotations
        .iter()
        .zip(operations.translations.iter())
    {
        let mut visited = vec![false; cell.num_atoms()];
        for i in 0..cell.num_atoms() {
            let new_pos = rotation.map(|e| e as f64) * cell.positions[i] + translation;
            let mut overlap = false;
            for j in 0..cell.num_atoms() {
                if visited[j] {
                    continue;
                }
                let mut diff = new_pos - cell.positions[j];
                diff -= diff.map(|x| x.round());
                if diff.iter().all(|x| x.abs() < 1e-4) {
                    visited[j] = true;
                    overlap = true;
                    break;
                }
            }
            if !overlap {
                return false;
            }
        }
    }

    true
}

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

    assert!(check_operations(&cell, &dataset.operations));

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

    assert!(check_operations(&cell, &dataset.operations));

    assert_eq!(dataset.number, 136); // P4_2/mnm
    assert_eq!(dataset.hall_number, 419);
    assert_eq!(dataset.num_operations(), 16);
    assert_eq!(dataset.orbits, vec![0, 0, 2, 2, 2, 2]);
}
