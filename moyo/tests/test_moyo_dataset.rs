#[macro_use]
extern crate approx;

use nalgebra::{matrix, vector, Matrix3, Vector3};
use serde_json;
use std::fs;
use std::path::Path;
use test_log::test;

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
    let num_operations = dataset.operations.len();
    for i in 0..num_operations {
        for j in i + 1..num_operations {
            if dataset.operations[i].rotation != dataset.operations[j].rotation {
                continue;
            }
            let mut diff = dataset.operations[i].translation - dataset.operations[j].translation;
            diff -= diff.map(|x| x.round());
            assert!(diff.iter().any(|x| x.abs() > 1e-4));
        }
    }

    for operation in dataset.operations.iter() {
        // Check if operation induces permutation
        let permutation =
            permutation_from_operation(cell, &operation.rotation, &operation.translation).unwrap();

        // For pure translation, check if mapped to the same site
        if operation.rotation == Rotation::identity() {
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
    assert_relative_eq!(
        prim_std_linear_inv,
        prim_std_linear_inv.map(|e| e.round()),
        epsilon = 1e-8
    );

    // Check std_rotation_matrix and std_linear
    assert_relative_eq!(
        dataset.std_rotation_matrix * cell.lattice.basis * dataset.std_linear,
        dataset.std_cell.lattice.basis,
        epsilon = 1e-8
    );
    // Check std_rotation_matrix and prim_std_linear
    assert_relative_eq!(
        dataset.std_rotation_matrix * cell.lattice.basis * dataset.prim_std_linear,
        dataset.prim_std_cell.lattice.basis,
        epsilon = 1e-8
    );
    // TODO: std_origin_shift
    // TODO: prim_origin_shift

    assert_eq!(dataset.mapping_std_prim.len(), cell.num_atoms());

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
    assert_eq!(dataset.wyckoffs, vec!['a', 'a', 'a', 'a']);
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
    assert_eq!(dataset.wyckoffs, vec!['a', 'a', 'f', 'f', 'f', 'f']);
}

#[test]
fn test_with_hcp() {
    // hcp, P6_3/mmc (No. 194)
    // https://next-gen.materialsproject.org/materials/mp-153
    let a = 3.17;
    let c = 5.14;
    let lattice = Lattice::new(matrix![
        a, 0.0, 0.0;
        -a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0;
        0.0, 0.0, c;
    ]);
    let positions = vec![
        // 2c
        Vector3::new(1.0 / 3.0, 2.0 / 3.0, 1.0 / 4.0),
        Vector3::new(2.0 / 3.0, 1.0 / 3.0, 3.0 / 4.0),
    ];
    let numbers = vec![0, 0];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 194);
    assert_eq!(dataset.hall_number, 488);
    assert_eq!(dataset.num_operations(), 24);
    assert_eq!(dataset.orbits, vec![0, 0]);
    // 2c and 2d belong to the same Wyckoff set
    assert_eq!(dataset.wyckoffs[0], dataset.wyckoffs[1]);
    if dataset.wyckoffs[0] != 'c' && dataset.wyckoffs[0] != 'd' {
        panic!("Unexpected Wyckoff letter: {}", dataset.wyckoffs[0]);
    }
}

#[test]
fn test_with_wurtzite() {
    // P6_3mc (No. 186)
    // https://next-gen.materialsproject.org/materials/mp-560588
    let a = 3.81;
    let c = 6.24;
    let lattice = Lattice::new(matrix![
        a, 0.0, 0.0;
        -a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0;
        0.0, 0.0, c;
    ]);
    let z1_2b = 0.00014;
    let z2_2b = 0.37486;
    let positions = vec![
        // 2b
        Vector3::new(1.0 / 3.0, 2.0 / 3.0, z1_2b),
        Vector3::new(2.0 / 3.0, 1.0 / 3.0, z1_2b + 0.5),
        // 2b
        Vector3::new(1.0 / 3.0, 2.0 / 3.0, z2_2b),
        Vector3::new(2.0 / 3.0, 1.0 / 3.0, z2_2b + 0.5),
    ];
    let numbers = vec![0, 0, 1, 1];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 186);
    assert_eq!(dataset.hall_number, 480);
    assert_eq!(dataset.num_operations(), 12);
    assert_eq!(dataset.orbits, vec![0, 0, 2, 2]);
    // 2a and 2b belong to the same Wyckoff set
    assert_eq!(dataset.wyckoffs[0], dataset.wyckoffs[1]);
    if dataset.wyckoffs[0] != 'a' && dataset.wyckoffs[0] != 'b' {
        panic!("Unexpected Wyckoff letter: {}", dataset.wyckoffs[0]);
    }
}

#[test]
fn test_with_corundum() {
    // Corundum structure, R-3c (No. 167)
    // https://materialsproject.org/materials/mp-1143
    let a = 4.80502783;
    let c = 13.11625361;
    let lattice = Lattice::new(matrix![
        a, 0.0, 0.0;
        -a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0;
        0.0, 0.0, c;
    ]);
    let positions = vec![
        // Al (12c)
        Vector3::new(0.33333333, 0.66666667, 0.81457067),
        Vector3::new(0.66666667, 0.33333333, 0.68542933),
        Vector3::new(0.00000000, 0.00000000, 0.64790400),
        Vector3::new(0.33333333, 0.66666667, 0.51876267),
        Vector3::new(0.00000000, 0.00000000, 0.14790400),
        Vector3::new(0.33333333, 0.66666667, 0.01876267),
        Vector3::new(0.66666667, 0.33333333, 0.98123733),
        Vector3::new(0.00000000, 0.00000000, 0.85209600),
        Vector3::new(0.66666667, 0.33333333, 0.48123733),
        Vector3::new(0.00000000, 0.00000000, 0.35209600),
        Vector3::new(0.33333333, 0.66666667, 0.31457067),
        Vector3::new(0.66666667, 0.33333333, 0.18542933),
        // O (18e)
        Vector3::new(0.36052117, 0.33333333, 0.58333333),
        Vector3::new(0.69385450, 0.69385450, 0.75000000),
        Vector3::new(0.97281217, 0.63947883, 0.58333333),
        Vector3::new(0.66666667, 0.02718783, 0.58333333),
        Vector3::new(0.00000000, 0.30614550, 0.75000000),
        Vector3::new(0.30614550, 0.00000000, 0.75000000),
        Vector3::new(0.02718783, 0.66666667, 0.91666667),
        Vector3::new(0.36052117, 0.02718783, 0.08333333),
        Vector3::new(0.63947883, 0.97281217, 0.91666667),
        Vector3::new(0.33333333, 0.36052117, 0.91666667),
        Vector3::new(0.66666667, 0.63947883, 0.08333333),
        Vector3::new(0.97281217, 0.33333333, 0.08333333),
        Vector3::new(0.69385450, 0.00000000, 0.25000000),
        Vector3::new(0.02718783, 0.36052117, 0.41666667),
        Vector3::new(0.30614550, 0.30614550, 0.25000000),
        Vector3::new(0.00000000, 0.69385450, 0.25000000),
        Vector3::new(0.33333333, 0.97281217, 0.41666667),
        Vector3::new(0.63947883, 0.66666667, 0.41666667),
    ];
    let numbers = vec![
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    ];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 167);
    assert_eq!(dataset.hall_number, 460); // Hexagonal setting
    assert_eq!(dataset.num_operations(), 36);
    assert_eq!(
        dataset.orbits,
        vec![
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12,
        ]
    );
    assert_eq!(
        dataset.wyckoffs,
        vec![
            'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'e', 'e', 'e', 'e', 'e',
            'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e', 'e',
        ]
    );
}

#[test]
#[allow(non_snake_case)]
fn test_with_hexagonal_Sc() {
    // P6_122 (No. 178)
    // https://legacy.materialsproject.org/materials/mp-601273/
    let a = 3.234;
    let c = 16.386;
    let lattice = Lattice::new(matrix![
        a, 0.0, 0.0;
        -a / 2.0, a * 3.0_f64.sqrt() / 2.0, 0.0;
        0.0, 0.0, c;
    ]);
    let x_6a = 0.4702;
    let positions = vec![
        // 6a
        Vector3::new(x_6a, 0.0, 0.0),
        Vector3::new(0.0, x_6a, 1.0 / 3.0),
        Vector3::new(-x_6a, -x_6a, 2.0 / 3.0),
        Vector3::new(-x_6a, 0.0, 0.5),
        Vector3::new(0.0, -x_6a, 5.0 / 6.0),
        Vector3::new(x_6a, x_6a, 1.0 / 6.0),
    ];
    let numbers = vec![0, 0, 0, 0, 0, 0];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 178);
    assert_eq!(dataset.hall_number, 472);
    assert_eq!(dataset.num_operations(), 12);
    assert_eq!(dataset.orbits, vec![0, 0, 0, 0, 0, 0]);
    assert_eq!(dataset.wyckoffs, vec!['a', 'a', 'a', 'a', 'a', 'a']);
}

#[test]
#[allow(non_snake_case)]
fn test_with_trigonal_Sc() {
    // https://next-gen.materialsproject.org/materials/mp-1055932
    let lattice = Lattice::new(matrix![
        -0.882444, 0.564392, -3.041088;
        -1.66822, -2.81974, -0.089223;
        -1.521212, 2.804144, -0.808507;
    ]);
    let positions = vec![Vector3::new(0.999917, 3.2999999999999996e-05, 1.7e-05)];
    let numbers = vec![0];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 1e-1; // This structure is distorted
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 166);
    assert_eq!(dataset.hall_number, 458);
    assert_eq!(dataset.num_operations(), 12); // Rhombohedral setting
    assert_eq!(dataset.orbits, vec![0]);
    if dataset.wyckoffs[0] != 'a' && dataset.wyckoffs[0] != 'b' {
        panic!("Unexpected Wyckoff letter: {}", dataset.wyckoffs[0]);
    }
}

#[test]
#[allow(non_snake_case)]
fn test_with_clathrate_Si() {
    // https://next-gen.materialsproject.org/materials/mp-1201492/
    // Pa-3
    let path = Path::new("tests/assets/mp-1201492.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 205);
    assert_eq!(dataset.hall_number, 501);
    assert_eq!(dataset.num_operations(), 24);
}

#[test]
fn test_with_mp_1197586() {
    // https://next-gen.materialsproject.org/materials/mp-1197586
    let path = Path::new("tests/assets/mp-1197586.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-3; // 1e-4 gives C2/m
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 194); // P6_3/mmc
    assert_eq!(dataset.hall_number, 488);
    assert_eq!(dataset.num_operations(), 24);
}

#[test]
fn test_with_mp_1185639() {
    // https://next-gen.materialsproject.org/materials/mp-1185639
    let path = Path::new("tests/assets/mp-1185639.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-2;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 187); // P-6m2
    assert_eq!(dataset.hall_number, 481);
    assert_eq!(dataset.num_operations(), 12);
}

#[test]
fn test_with_mp_1221598() {
    // https://next-gen.materialsproject.org/materials/mp-1221598
    let path = Path::new("tests/assets/mp-1221598.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-1;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 225); // Fm-3m
}

#[test]
fn test_with_mp_569901() {
    let path = Path::new("tests/assets/mp-569901.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-1;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);

    assert_eq!(dataset.number, 118); // P-4n2
}

#[test]
fn test_with_mp_30665() {
    let path = Path::new("tests/assets/mp-30665.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-1;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);
}

#[test]
fn test_monotonic_symmetry_recovery() {
    let path = Path::new("tests/assets/mp-1277787.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let symprec = 5e-2;
    let dataset1 = MoyoDataset::new(&cell, symprec, angle_tolerance, setting).unwrap();

    let symprec = 1e-1;
    let dataset2 = MoyoDataset::new(&cell, symprec, angle_tolerance, setting).unwrap();

    assert!(dataset1.number <= dataset2.number);
}

#[test]
fn test_with_mp_550745() {
    let path = Path::new("tests/assets/mp-550745.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-2;
    let angle_tolerance = AngleTolerance::Radian(0.1);
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);
}

#[test]
fn test_niggli_reduction_corner_cases() {
    // https://github.com/spglib/moyo/issues/35
    let symprec = 1e-5;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    for path in vec![
        Path::new("tests/assets/wbm-1-42389.json"),
        Path::new("tests/assets/wbm-1-42433.json"),
    ] {
        let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

        let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
        assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
        assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);
    }
}

#[test]
fn test_primitive_symmetry_search_corner_case() {
    // https://github.com/spglib/moyo/issues/38
    let path = Path::new("tests/assets/wbm-1-29497.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-2;
    let angle_tolerance = AngleTolerance::Default;
    let setting = Setting::Standard;

    let dataset = assert_dataset(&cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.std_cell, symprec, angle_tolerance, setting);
    assert_dataset(&dataset.prim_std_cell, symprec, angle_tolerance, setting);
}

#[test]
fn test_with_high_symprec_and_angle_tolerance() {
    let lattice = Lattice::new(Matrix3::identity());
    let positions = vec![
        vector![0.0, 0.0, 0.0],
        vector![0.0, 0.5, 0.5],
        vector![0.5, 0.0, 0.5],
        vector![0.5, 0.5, 0.0],
    ];
    let numbers = vec![0, 0, 0, 0];
    let cell = Cell::new(lattice, positions, numbers);

    let symprec = 0.1;
    let angle_tolerance = AngleTolerance::Radian(1.0);
    let setting = Setting::Spglib;

    let _ = MoyoDataset::new(&cell, symprec, angle_tolerance, setting).unwrap();
}
