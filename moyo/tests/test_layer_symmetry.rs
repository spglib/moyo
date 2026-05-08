use nalgebra::{Vector3, matrix};
use test_log::test;

use moyo::base::{AngleTolerance, Cell, Lattice, LayerCell, Rotation};
use moyo::search::LayerPrimitiveSymmetrySearch;

const SYMPREC: f64 = 1e-4;

/// Verify the layer-group block form (paper eq. 4):
///   W[0,2] = W[1,2] = W[2,0] = W[2,1] = 0 and |W[2,2]| = 1.
fn assert_layer_block_form(rotation: &Rotation) {
    assert_eq!(rotation[(0, 2)], 0);
    assert_eq!(rotation[(1, 2)], 0);
    assert_eq!(rotation[(2, 0)], 0);
    assert_eq!(rotation[(2, 1)], 0);
    assert_eq!(rotation[(2, 2)].abs(), 1);
}

/// LG 2 (p-1): triclinic-oblique in-plane lattice with two atoms accidentally
/// placed at inversion-related positions through (0.2, 0.35, 0.15). The
/// inversion has a non-zero t_z compensating for its off-origin center.
#[test]
fn test_layer_p_minus_1_round_trip() {
    // a = 1.0 along x, b = 1.5 at gamma = 80 deg, c = 5.0 along z.
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

    let layer = LayerCell::new(cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();
    let result =
        LayerPrimitiveSymmetrySearch::new(&layer, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    // Order of LG 2 is 2 (identity + inversion).
    assert_eq!(result.operations.len(), 2);
    let inversion: Rotation = matrix![-1, 0, 0; 0, -1, 0; 0, 0, -1];
    let has_identity = result
        .operations
        .iter()
        .any(|op| op.rotation == Rotation::identity());
    let has_inversion = result.operations.iter().any(|op| op.rotation == inversion);
    assert!(has_identity, "LG p-1 must contain the identity");
    assert!(has_inversion, "LG p-1 must contain the inversion");
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }
}

/// LG 6 (p2/m): monoclinic-oblique with both atoms on z = 0.07. The shared
/// z plane carries a horizontal mirror (m_z), which combined with the
/// 2-fold along c lifts the group from p112 to p2/m.
#[test]
fn test_layer_p2_per_m_round_trip() {
    let gamma = 80.0_f64.to_radians();
    let a = 1.0;
    let b = 1.5;
    let cell = Cell::new(
        Lattice::new(matrix![
            a, b * gamma.cos(), 0.0;
            0.0, b * gamma.sin(), 0.0;
            0.0, 0.0, 5.0;
        ]),
        vec![
            Vector3::new(0.31, 0.42, 0.07),
            Vector3::new(-0.31, -0.42, 0.07),
        ],
        vec![1, 1],
    );

    let layer = LayerCell::new(cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();
    let result =
        LayerPrimitiveSymmetrySearch::new(&layer, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    // Order of LG 6 (p2/m) is 4: {1, 2_c, m_z, -1}.
    assert_eq!(result.operations.len(), 4);
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }

    // Spot-check the generators: 2-fold along c and the horizontal mirror m_z.
    let two_fold_c: Rotation = matrix![-1, 0, 0; 0, -1, 0; 0, 0, 1];
    let mirror_z: Rotation = matrix![1, 0, 0; 0, 1, 0; 0, 0, -1];
    let has_two_fold = result.operations.iter().any(|op| op.rotation == two_fold_c);
    let has_mirror_z = result.operations.iter().any(|op| op.rotation == mirror_z);
    assert!(has_two_fold, "LG p2/m must contain the 2-fold along c");
    assert!(
        has_mirror_z,
        "LG p2/m must contain the horizontal mirror m_z"
    );
}

/// LG 61 (p4/mmm): square in-plane lattice with a single atom at (0,0,0.1).
/// The atom's own z plane carries a horizontal mirror, so the layer group
/// is the centrosymmetric supergroup p4/mmm of order 16, not p4mm.
#[test]
fn test_layer_p4_per_mmm_round_trip() {
    let cell = Cell::new(
        Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]),
        vec![Vector3::new(0.0, 0.0, 0.1)],
        vec![1],
    );

    let layer = LayerCell::new(cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();
    let result =
        LayerPrimitiveSymmetrySearch::new(&layer, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    assert_eq!(result.operations.len(), 16);
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }

    // Spot-check generators: 4-fold along c, mirror across x = 0, and m_z.
    let four_fold_c: Rotation = matrix![0, -1, 0; 1, 0, 0; 0, 0, 1];
    let mirror_x: Rotation = matrix![-1, 0, 0; 0, 1, 0; 0, 0, 1];
    let mirror_z: Rotation = matrix![1, 0, 0; 0, 1, 0; 0, 0, -1];
    assert!(
        result
            .operations
            .iter()
            .any(|op| op.rotation == four_fold_c),
        "LG p4/mmm must contain the 4-fold along c"
    );
    assert!(
        result.operations.iter().any(|op| op.rotation == mirror_x),
        "LG p4/mmm must contain the mirror across x = 0"
    );
    assert!(
        result.operations.iter().any(|op| op.rotation == mirror_z),
        "LG p4/mmm must contain the horizontal mirror m_z"
    );
}
