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

/// RhO2 monolayer (JARVIS jid=JVASP-77581), bulk parent space group P-3m1
/// (no. 164). The layer is centrosymmetric (D_3d): inversion through the
/// origin maps Rh -> Rh and the two O atoms onto each other. moyo
/// correctly recovers the full LG 72 (p-3m1) of order 12; spglib's
/// `get_symmetry_layerdataset` returns only the c-fixing C_3v subset for
/// this structure (number=69, 6 ops), missing the inversion. This test
/// pins moyo to the centrosymmetric answer.
#[test]
fn test_layer_dataset_rho2_jvasp_77581() {
    let cell = Cell::new(
        Lattice::new(matrix![
            3.077561866713373,   -5.5895361638e-06,  0.0;
            -1.538785277570888,  2.665244603795193,  0.0;
            0.0,                 0.0,                21.887078;
        ]),
        vec![
            Vector3::new(3.98962103e-08, -1.0507825500000001e-08, 2.690075935e-07),
            Vector3::new(0.33333345542078574, 0.6666665498487259, 0.0424523285784089),
            Vector3::new(0.6666665046829969, 0.33333346065909963, 0.9575474024139978),
        ],
        vec![45, 8, 8],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();
    assert_eq!(dataset.number, 72);
    assert_eq!(dataset.num_operations(), 12);
}

/// Repro for the moyopy "Wyckoff position assignment failed" crash on a
/// centrosymmetric LG 72 (p-3m1) layer with the inversion center placed at
/// z = 0.3 (rather than z = 0). The structure is the same as RhO2 above
/// modulo an origin shift of z by 0.3; the layer group is unchanged.
#[test]
fn test_layer_dataset_rho2_off_origin() {
    let a = 3.077561866713373;
    let cell = Cell::new(
        Lattice::new(matrix![
            a,        0.0,                              0.0;
            -0.5 * a, a * 0.866025403784438646763_f64,  0.0;
            0.0,      0.0,                              21.887078;
        ]),
        vec![
            Vector3::new(0.0, 0.0, 0.3),
            Vector3::new(2.0 / 3.0, 1.0 / 3.0, 0.3 + 0.0425),
            Vector3::new(1.0 / 3.0, 2.0 / 3.0, 0.3 - 0.0425),
        ],
        vec![45, 8, 8],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();
    assert_eq!(dataset.number, 72);
}

/// Bi2Pt monolayer (JARVIS jid=JVASP-76516), bulk parent SG 164 P-3m1.
/// The inversion center sits at the Pt site, z = 0.0626 -- not at z = 0
/// nor at z = 1/2. spglib finds LG 72 / 12 ops; bulk MoyoDataset finds
/// SG 164 / 12 ops; the layer identifier should agree.
#[test]
fn test_layer_dataset_bi2pt_jvasp_76516() {
    let cell = Cell::new(
        Lattice::new(matrix![
            3.631374852245063,    -2.0965748614627038, 0.0;
            0.0,                   4.1931497229254076, 0.0;
            0.0,                   0.0,                22.860029;
        ]),
        vec![
            Vector3::new(0.33333300000000315, 0.6666669999999968, 0.122798185664739),
            Vector3::new(0.6666669999999968, 0.3333330000000032, 0.0023118143352603),
            Vector3::new(0.0, 0.0, 0.0625550000000032),
        ],
        vec![83, 83, 78],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();
    assert_eq!(dataset.number, 72);
    assert_eq!(dataset.num_operations(), 12);
}

/// SiP monolayer (JARVIS jid=JVASP-27741), bulk parent SG 12 C2/m.
/// spglib finds LG 18 / 8 ops; moyo's `LayerGroup::new` currently fails
/// because the C-centered primitive-basis ambiguity needs a
/// centering-aware correction (the GL(2, Z) shears between equivalent
/// primitive cells of an oc Bravais lattice are infinite-order, so a
/// static finite list of correction matrices cannot cover them).
/// Marked `#[ignore]` until the centering-aware correction
/// machinery (analogue of bulk `correction_transformation_matrices`)
/// lands. JARVIS bench bucket: dominant in SG 12 / 21 / 65.
#[test]
#[ignore]
fn test_layer_dataset_sip_jvasp_27741() {
    let cell = Cell::new(
        Lattice::new(matrix![
            3.544338948503597,    0.0,                  0.0;
            0.0,                  20.623454619019398,   0.0;
            0.0,                  0.0,                  27.611119;
        ]),
        vec![
            Vector3::new(0.5, 0.7018568074325104, 0.099345553246816),
            Vector3::new(0.5, 0.0670975874506597, 0.193051629621707),
            Vector3::new(0.5, 0.7051327661543786, 0.18505817140897843),
            Vector3::new(0.0, 0.20513276615437862, 0.18505817140897843),
            Vector3::new(0.0, 0.9397388688497569, 0.1571337073357672),
            Vector3::new(0.5, 0.43973886884975705, 0.1571337073357672),
            Vector3::new(0.0, 0.5670975874506596, 0.193051629621707),
            Vector3::new(0.0, 0.8292144252409636, 0.1352653645591983),
            Vector3::new(0.5, 0.0638204858451885, 0.1073385088143841),
            Vector3::new(0.0, 0.5638204858451881, 0.1073385088143841),
            Vector3::new(0.0, 0.20185680743251064, 0.099345553246816),
            Vector3::new(0.5, 0.32921442524096395, 0.1352653645591983),
            Vector3::new(0.0, 0.4582284213727509, 0.208414460112573),
            Vector3::new(0.5, 0.9582284213727511, 0.208414460112573),
            Vector3::new(0.5, 0.26954783797573045, 0.2064242218219204),
            Vector3::new(0.0, 0.7695478379757302, 0.2064242218219204),
            Vector3::new(0.5, 0.49940495566583837, 0.0859734338638388),
            Vector3::new(0.0, 0.3107261843718542, 0.08398562506008961),
            Vector3::new(0.5, 0.8107261843718538, 0.08398562506008961),
            Vector3::new(0.0, 0.6549806360416701, 0.0603195954060404),
            Vector3::new(0.5, 0.15498063604167003, 0.0603195954060404),
            Vector3::new(0.0, 0.1139740235986913, 0.23207872874868993),
            Vector3::new(0.0, 0.9994049556658384, 0.0859734338638388),
            Vector3::new(0.5, 0.6139740235986916, 0.23207872874868993),
        ],
        vec![
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
            15, 15,
        ],
    );

    let dataset = MoyoLayerDataset::with_default(&cell, SYMPREC).unwrap();
    assert_eq!(dataset.number, 18);
    assert_eq!(dataset.num_operations(), 8);
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
