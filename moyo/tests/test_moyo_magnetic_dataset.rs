#[macro_use]
extern crate approx;

use nalgebra::{matrix, vector, Matrix3};
use test_log::test;

use moyo::base::{
    AngleTolerance, Collinear, Lattice, MagneticCell, MagneticMoment, RotationMagneticMomentAction,
};
use moyo::MoyoMagneticDataset;

/// Sanity-check MoyoMagneticDataset
fn assert_magnetic_dataset<M: MagneticMoment>(
    magnetic_cell: &MagneticCell<M>,
    symprec: f64,
    angle_tolerance: AngleTolerance,
    mag_symprec: Option<f64>,
    action: RotationMagneticMomentAction,
) -> MoyoMagneticDataset<M> {
    let dataset =
        MoyoMagneticDataset::new(magnetic_cell, symprec, angle_tolerance, mag_symprec, action)
            .unwrap();

    // std_mag_cell
    let std_dataset = MoyoMagneticDataset::new(
        &dataset.std_mag_cell,
        symprec,
        angle_tolerance,
        mag_symprec,
        action,
    )
    .unwrap();
    assert_eq!(std_dataset.uni_number, dataset.uni_number);

    // prim_std_mag_cell
    let prim_std_dataset = MoyoMagneticDataset::new(
        &dataset.prim_std_mag_cell,
        symprec,
        angle_tolerance,
        mag_symprec,
        action,
    )
    .unwrap();
    assert_eq!(prim_std_dataset.uni_number, dataset.uni_number);

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
        dataset.std_rotation_matrix * magnetic_cell.cell.lattice.basis * dataset.std_linear,
        dataset.std_mag_cell.cell.lattice.basis,
        epsilon = 1e-8
    );
    // Check std_rotation_matrix and prim_std_linear
    assert_relative_eq!(
        dataset.std_rotation_matrix * magnetic_cell.cell.lattice.basis * dataset.prim_std_linear,
        dataset.prim_std_mag_cell.cell.lattice.basis,
        epsilon = 1e-8
    );
    // TODO: std_origin_shift
    // TODO: prim_origin_shift

    assert_eq!(dataset.mapping_std_prim.len(), magnetic_cell.num_atoms());

    dataset
}

#[test]
fn test_with_rutile() {
    let lattice = Lattice::new(Matrix3::identity());
    let positions = vec![
        // Ti (2a)
        vector![0.0, 0.0, 0.0],
        vector![0.5, 0.5, 0.5],
        // O (4f)
        vector![0.3, 0.3, 0.0],
        vector![0.7, 0.7, 0.0],
        vector![0.2, 0.8, 0.5],
        vector![0.8, 0.2, 0.5],
    ];
    let numbers = vec![0, 0, 1, 1, 1, 1];

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let mag_symprec = None;
    let action = RotationMagneticMomentAction::Polar;

    {
        // Type-I, 136.495: -P 4n 2n
        let magmoms = vec![
            Collinear(0.7),
            Collinear(0.7),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
        ];
        let magnetic_cell =
            MagneticCell::new(lattice.clone(), positions.clone(), numbers.clone(), magmoms);
        let dataset = assert_magnetic_dataset(
            &magnetic_cell,
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
        );

        assert_eq!(dataset.uni_number, 1155);
    }

    {
        // Type-II, "136.496": -P 4n 2n 1'
        let magmoms = vec![
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
        ];
        let magnetic_cell =
            MagneticCell::new(lattice.clone(), positions.clone(), numbers.clone(), magmoms);
        let dataset = assert_magnetic_dataset(
            &magnetic_cell,
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
        );

        assert_eq!(dataset.uni_number, 1156);
    }

    {
        // Type-III, "136.498": -P 4n' 2n'
        let magmoms = vec![
            Collinear(0.7),
            Collinear(-0.7),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
        ];
        let magnetic_cell =
            MagneticCell::new(lattice.clone(), positions.clone(), numbers.clone(), magmoms);
        let dataset = assert_magnetic_dataset(
            &magnetic_cell,
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
        );

        assert_eq!(dataset.uni_number, 1158);
        assert_eq!(dataset.num_magnetic_operations(), 16);
        assert_eq!(dataset.orbits, vec![0, 0, 2, 2, 2, 2]);
        assert_eq!(dataset.std_mag_cell.num_atoms(), 6);
        assert_eq!(dataset.prim_std_mag_cell.num_atoms(), 6);
        assert_eq!(dataset.mapping_std_prim, vec![0, 1, 2, 3, 4, 5]);
    }
}

#[test]
fn test_with_rutile_type4() {
    let lattice = Lattice::new(matrix![
        5.0, 0.0, 0.0;
        0.0, 5.0, 0.0;
        0.0, 0.0, 6.0;
    ]);
    let positions = vec![
        // Ti (2a)
        vector![0.0, 0.0, 0.0],
        vector![0.5, 0.5, 0.25],
        // O (4f)
        vector![0.3, 0.3, 0.0],
        vector![0.7, 0.7, 0.0],
        vector![0.2, 0.8, 0.25],
        vector![0.8, 0.2, 0.25],
        // Ti (2a)
        vector![0.0, 0.0, 0.5],
        vector![0.5, 0.5, 0.75],
        // O (4f)
        vector![0.3, 0.3, 0.5],
        vector![0.7, 0.7, 0.5],
        vector![0.2, 0.8, 0.75],
        vector![0.8, 0.2, 0.75],
    ];
    let numbers = vec![0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1];
    let magmoms = vec![
        // Ti (2a)
        Collinear(0.3),
        Collinear(0.3),
        // O (4f)
        Collinear(0.0),
        Collinear(0.0),
        Collinear(0.0),
        Collinear(0.0),
        // Ti (2a)
        Collinear(-0.3),
        Collinear(-0.3),
        // O (4f)
        Collinear(0.0),
        Collinear(0.0),
        Collinear(0.0),
        Collinear(0.0),
    ];
    let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magmoms);

    let symprec = 1e-4;
    let angle_tolerance = AngleTolerance::Default;
    let mag_symprec = None;
    let action = RotationMagneticMomentAction::Polar;

    let dataset = assert_magnetic_dataset(
        &magnetic_cell,
        symprec,
        angle_tolerance,
        mag_symprec,
        action,
    );

    assert_eq!(dataset.uni_number, 932);
}
