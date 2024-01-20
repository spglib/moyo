#[macro_use]
extern crate approx;

use nalgebra::{matrix, vector, Matrix3, Vector3};

use moyo::base::{AngleTolerance, Cell, Lattice, Permutation, Rotation, Translation};
use moyo::data::Setting;
use moyo::MoyoDataset;

/// Sanity-check MoyoDataset
fn assert_dataset(
    cell: &Cell,
    symprec: f64,
    angle_tolerance: AngleTolerance,
    setting: Setting,
) -> MoyoDataset {
    let dataset = MoyoDataset::new(cell, symprec, angle_tolerance, setting).unwrap();

    // Check if operations are unique
    let num_operations = dataset.operations.num_operations();
    for i in 0..num_operations {
        for j in i + 1..num_operations {
            if dataset.operations.rotations[i] != dataset.operations.rotations[j] {
                continue;
            }
            let mut diff = dataset.operations.translations[i] - dataset.operations.translations[j];
            diff -= diff.map(|x| x.round());
            assert!(diff.iter().any(|x| x.abs() > 1e-4));
        }
    }

    for (rotation, translation) in dataset
        .operations
        .rotations
        .iter()
        .zip(dataset.operations.translations.iter())
    {
        // Check if operation induces permutation
        let permutation = permutation_from_operation(cell, rotation, translation).unwrap();

        // For pure translation, check if mapped to the same site
        if rotation == &Rotation::identity() {
            for i in 0..cell.num_atoms() {
                let j = permutation.apply(i);
                assert_eq!(dataset.mapping_std_prim[i], dataset.mapping_std_prim[j]);
            }
        }

        for i in 0..cell.num_atoms() {
            let j = permutation.apply(i);
            assert_eq!(cell.numbers[i], cell.numbers[j]);
            // Check if belong to the same orbit
            assert_eq!(dataset.orbits[i], dataset.orbits[j]);
        }
    }

    // std_cell
    let std_dataset =
        MoyoDataset::new(&dataset.std_cell, symprec, angle_tolerance, setting).unwrap();
    assert_eq!(std_dataset.number, dataset.number);
    assert_eq!(std_dataset.hall_number, dataset.hall_number);

    // prim_std_cell
    let prim_std_dataset =
        MoyoDataset::new(&dataset.prim_std_cell, symprec, angle_tolerance, setting).unwrap();
    assert_eq!(prim_std_dataset.number, dataset.number);
    assert_eq!(prim_std_dataset.hall_number, dataset.hall_number);

    // prim_std_linear should be an inverse of an integer matrix
    let prim_std_linear_inv = dataset
        .prim_std_linear
        .map(|e| e as f64)
        .try_inverse()
        .unwrap();
    assert_relative_eq!(prim_std_linear_inv, prim_std_linear_inv.map(|e| e.round()));

    // Check std_rotation_matrix and std_linear
    assert_relative_eq!(
        dataset.std_rotation_matrix * cell.lattice.basis * dataset.std_linear,
        dataset.std_cell.lattice.basis
    );
    // Check std_rotation_matrix and prim_std_linear
    assert_relative_eq!(
        dataset.std_rotation_matrix * cell.lattice.basis * dataset.prim_std_linear,
        dataset.prim_std_cell.lattice.basis
    );
    // TODO: std_origin_shift
    // TODO: prim_origin_shift

    dataset
}

/// O(num_atoms^2)
fn permutation_from_operation(
    cell: &Cell,
    rotation: &Rotation,
    translation: &Translation,
) -> Option<Permutation> {
    let mut visited = vec![false; cell.num_atoms()];
    let mut mapping = vec![0; cell.num_atoms()];
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
                mapping[i] = j;
                overlap = true;
                break;
            }
        }
        if !overlap {
            return None;
        }
    }
    Some(Permutation::new(mapping))
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
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Spglib;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

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
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Spglib;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 136); // P4_2/mnm
    assert_eq!(dataset.hall_number, 419);
    assert_eq!(dataset.num_operations(), 16);
    assert_eq!(dataset.orbits, vec![0, 0, 2, 2, 2, 2]);
}
