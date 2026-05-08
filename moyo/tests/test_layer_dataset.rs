#[macro_use]
extern crate approx;

use nalgebra::{Vector3, matrix};
use test_log::test;

use moyo::MoyoLayerDataset;
use moyo::base::{AngleTolerance, Cell, Lattice};
use moyo::data::LayerSetting;

const SYMPREC: f64 = 1e-4;

/// LG 2 (p-1): triclinic-oblique with two atoms accidentally placed at
/// inversion-related positions through (0.2, 0.35, 0.15). Verifies
/// identification and standardized-cell convention end-to-end through the
/// public dataset.
#[test]
fn test_layer_dataset_p_minus_1() {
    let gamma = 80.0_f64.to_radians();
    let a = 1.0;
    let b = 1.5;
    let cell = Cell::new(
        Lattice::new(matrix![
            a, b * gamma.cos(), 0.0;
            0.0, b * gamma.sin(), 0.0;
            0.0, 0.0, 5.0;
        ]),
        vec![Vector3::new(0.1, 0.2, 0.1), Vector3::new(0.3, 0.5, 0.2)],
        vec![1, 1],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();

    assert_eq!(dataset.number, 2);
    assert_eq!(dataset.hall_number, 2);
    // Order of LG 2 is 2; the inversion has a non-zero t_z compensating for
    // the off-origin inversion center.
    assert_eq!(dataset.num_operations(), 2);
    assert_eq!(dataset.orbits, vec![0, 0]);
    assert_eq!(dataset.wyckoffs, vec!['e', 'e']);
    // Standardized output: |c| preserved, c along z, a along x.
    let std_basis = dataset.std_cell.lattice().basis();
    assert_relative_eq!(std_basis.column(2).norm(), 5.0, epsilon = 1e-10);
    assert_relative_eq!(std_basis[(2, 2)], 5.0, epsilon = 1e-10);
    assert_relative_eq!(std_basis[(0, 2)], 0.0, epsilon = 1e-10);
    assert_relative_eq!(std_basis[(1, 2)], 0.0, epsilon = 1e-10);
    // mp Bravais + 2 atoms.
    assert_eq!(dataset.pearson_symbol, "mp2");
}

/// LG 6 (p2/m): primitive monoclinic-oblique, two atoms placed at the same
/// z plane and related by the 2-fold along c. The shared z=0.07 plane also
/// hosts a horizontal mirror, lifting the group from p112 to p2/m. Verifies
/// orbit collapse and wyckoff equality across the centro-symmetric pair.
#[test]
fn test_layer_dataset_p2_per_m_orbit_collapse() {
    let cell = Cell::new(
        Lattice::new(matrix![
            1.0, 0.2, 0.0;
            0.0, 1.3, 0.0;
            0.0, 0.0, 5.0;
        ]),
        vec![
            Vector3::new(0.31, 0.42, 0.07),
            Vector3::new(-0.31, -0.42, 0.07),
        ],
        vec![1, 1],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();

    assert_eq!(dataset.number, 6);
    // Two atoms related by the 2-fold along c -> one orbit.
    assert_eq!(dataset.orbits, vec![0, 0]);
    assert_eq!(dataset.wyckoffs[0], dataset.wyckoffs[1]);
    // mp + 2 atoms in the conventional cell.
    assert_eq!(dataset.pearson_symbol, "mp2");
}

/// LG 61 (p4/mmm): square in-plane lattice with a single atom at (0,0,z).
/// The atom's own z plane carries a horizontal mirror, so the layer group
/// is p4/mmm (the centro-symmetric supergroup of p4mm) rather than LG 55.
/// Verifies std_cell snaps to a square (a = b) and the site symmetry is 4/mmm.
#[test]
fn test_layer_dataset_p4_per_mmm() {
    let cell = Cell::new(
        Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]),
        vec![Vector3::new(0.0, 0.0, 0.1)],
        vec![1],
    );

    let dataset = MoyoLayerDataset::new(
        &cell,
        SYMPREC,
        AngleTolerance::default(),
        LayerSetting::Standard,
        true,
    )
    .unwrap();

    assert_eq!(dataset.number, 61);
    assert_eq!(dataset.num_operations(), 16);
    // Standardized cell is square with c along z.
    let basis = dataset.std_cell.lattice().basis();
    let a = basis.column(0).norm();
    let b = basis.column(1).norm();
    assert_relative_eq!(a, b, epsilon = 1e-10);
    assert_relative_eq!(basis.column(2).norm(), 5.0, epsilon = 1e-10);
    // Single high-symmetry orbit (Wyckoff 'a', site symmetry is 4/mmm).
    assert_eq!(dataset.wyckoffs, vec!['a']);
    assert!(dataset.site_symmetry_symbols[0].contains("4/mmm"));
    // tp Bravais + 1 atom.
    assert_eq!(dataset.pearson_symbol, "tp1");
}

/// CoO2 monolayer (JARVIS jid=JVASP-14441), exfoliated 2D structure with
/// c orthogonal to the in-plane vectors. spglib reports LG 72; this test
/// pins moyo to the same answer.
#[test]
fn test_layer_dataset_coo2_jvasp_14441() {
    let cell = Cell::new(
        Lattice::new(matrix![
            2.8173501021142346,  0.0,                0.0;
            -1.4086750510571173, 2.439896157530945,  0.0;
            0.0,                 0.0,                24.551775;
        ]),
        vec![
            Vector3::new(0.2690539999999971, 0.9290120000000003, 0.0926970000000011),
            Vector3::new(0.6023870000000002, 0.595678999999997, 0.1307621356947131),
            Vector3::new(0.9357210000000009, 0.26234500000000344, 0.0546328643052866),
        ],
        vec![27, 8, 8],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();
    assert_eq!(dataset.number, 72);
}

/// Tilted-c input must be rejected by the layer-cell constructor.
#[test]
fn test_layer_dataset_rejects_tilted_c() {
    let cell = Cell::new(
        Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.5, 0.0, 5.0;
        ]),
        vec![Vector3::new(0.0, 0.0, 0.0)],
        vec![1],
    );
    let result = MoyoLayerDataset::with_default(&cell, SYMPREC);
    assert!(result.is_err());
}
