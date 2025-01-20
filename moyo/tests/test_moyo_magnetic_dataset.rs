use nalgebra::{vector, Matrix3};
use test_log::test;

use moyo::base::{AngleTolerance, Collinear, Lattice, MagneticCell, RotationMagneticMomentAction};
use moyo::MoyoMagneticDataset;

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
        // Type-III, "136.498": -P 4n' 2n'
        let magmoms = vec![
            Collinear(0.7),
            Collinear(-0.7),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
            Collinear(0.0),
        ];
        let magnetic_cell = MagneticCell::new(lattice, positions, numbers, magmoms);
        let dataset = MoyoMagneticDataset::new(
            &magnetic_cell,
            symprec,
            angle_tolerance,
            mag_symprec,
            action,
        )
        .unwrap();

        assert_eq!(dataset.uni_number, 1158);
    }
}
