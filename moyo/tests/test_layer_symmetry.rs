use nalgebra::{Vector3, matrix};
use test_log::test;

use moyo::base::{AngleTolerance, Cell, Lattice, Rotation};
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

/// LG 1 (p1): triclinic-oblique in-plane lattice, atoms placed so that no
/// non-trivial Bravais element survives. Only the identity is a symmetry.
#[test]
fn test_layer_p1_round_trip() {
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
        // WHY: two atoms at unrelated generic positions kill 2_c, -1, and m_z
        // simultaneously -- a single atom would always admit the full holohedry
        // because translations are free for one site.
        vec![Vector3::new(0.1, 0.2, 0.1), Vector3::new(0.3, 0.5, 0.2)],
        vec![1, 1],
    );

    let result =
        LayerPrimitiveSymmetrySearch::new(&cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    // Order of LG 1 is 1 (identity only).
    assert_eq!(result.operations.len(), 1);
    assert_eq!(result.operations[0].rotation, Rotation::identity());
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }
}

/// LG 3 (p112): monoclinic-oblique with the 2-fold along c (the aperiodic axis).
/// Generic atom position picks up only identity + 2-fold along c.
#[test]
fn test_layer_p2_round_trip() {
    let gamma = 80.0_f64.to_radians();
    let a = 1.0;
    let b = 1.5;
    let cell = Cell::new(
        Lattice::new(matrix![
            a, b * gamma.cos(), 0.0;
            0.0, b * gamma.sin(), 0.0;
            0.0, 0.0, 5.0;
        ]),
        // WHY: generic z != 0 still respects the 2-fold along c (z -> z).
        vec![
            Vector3::new(0.31, 0.42, 0.07),
            Vector3::new(-0.31, -0.42, 0.07),
        ],
        vec![1, 1],
    );

    let result =
        LayerPrimitiveSymmetrySearch::new(&cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    // Order of LG 3 (p112) is 2.
    assert_eq!(result.operations.len(), 2);
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }

    // The 2-fold along c is diag(-1, -1, 1).
    let two_fold_c: Rotation = matrix![-1, 0, 0; 0, -1, 0; 0, 0, 1];
    let has_two_fold = result.operations.iter().any(|op| op.rotation == two_fold_c);
    assert!(has_two_fold, "LG p112 must contain the 2-fold along c");
}

/// LG 55 (p4mm): square in-plane lattice, atom at the 4mm site (origin) but
/// shifted to z != 0, ±1/2 to break m_z (which would promote 4mm to 4/mmm).
/// Order of 4mm is 8; all elements already fix the c-axis so the LG block-form
/// filter is a no-op here.
#[test]
fn test_layer_p4mm_round_trip() {
    let cell = Cell::new(
        Lattice::new(matrix![
            1.0, 0.0, 0.0;
            0.0, 1.0, 0.0;
            0.0, 0.0, 5.0;
        ]),
        // WHY: z = 0.1 is fixed by 4mm but not by m_z, isolating LG 55 from LG 61.
        vec![Vector3::new(0.0, 0.0, 0.1)],
        vec![1],
    );

    let result =
        LayerPrimitiveSymmetrySearch::new(&cell, SYMPREC, AngleTolerance::Radian(1e-2)).unwrap();

    assert_eq!(result.operations.len(), 8);
    for op in result.operations.iter() {
        assert_layer_block_form(&op.rotation);
    }

    // Spot-check generators: 4-fold along c and a mirror across x = 0.
    let four_fold_c: Rotation = matrix![0, -1, 0; 1, 0, 0; 0, 0, 1];
    let mirror_x: Rotation = matrix![-1, 0, 0; 0, 1, 0; 0, 0, 1];
    assert!(
        result
            .operations
            .iter()
            .any(|op| op.rotation == four_fold_c),
        "LG p4mm must contain the 4-fold along c"
    );
    assert!(
        result.operations.iter().any(|op| op.rotation == mirror_x),
        "LG p4mm must contain the mirror across x = 0"
    );
}
