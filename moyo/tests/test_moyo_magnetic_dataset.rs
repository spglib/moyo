#[macro_use]
extern crate approx;

use nalgebra::{Matrix3, Rotation3, matrix, vector};
use std::fs;
use std::path::Path;
use test_log::test;

use moyo::MoyoMagneticDataset;
use moyo::base::{
    AngleTolerance, Collinear, Lattice, MagneticCell, MagneticMoment, NonCollinear,
    RotationMagneticMomentAction,
};
use moyo::data::magnetic_operations_from_uni_number;

fn assert_magnetic_dataset_with_default<M: MagneticMoment>(
    magnetic_cell: &MagneticCell<M>,
    symprec: f64,
    action: RotationMagneticMomentAction,
) -> MoyoMagneticDataset<M> {
    assert_magnetic_dataset(
        magnetic_cell,
        symprec,
        AngleTolerance::default(),
        None,
        action,
    )
}

/// Sanity-check MoyoMagneticDataset
fn assert_magnetic_dataset<M: MagneticMoment>(
    magnetic_cell: &MagneticCell<M>,
    symprec: f64,
    angle_tolerance: AngleTolerance,
    mag_symprec: Option<f64>,
    action: RotationMagneticMomentAction,
) -> MoyoMagneticDataset<M> {
    let rotate_basis = true;
    let dataset = MoyoMagneticDataset::new(
        magnetic_cell,
        symprec,
        angle_tolerance,
        mag_symprec,
        action,
        rotate_basis,
    )
    .unwrap();

    // std_mag_cell
    let std_dataset = MoyoMagneticDataset::new(
        &dataset.std_mag_cell,
        symprec,
        angle_tolerance,
        mag_symprec,
        action,
        rotate_basis,
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
        rotate_basis,
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

    // Lattice refinement is a symmetric positive stretch before rotation.
    let stretch = dataset.std_rotation_matrix.transpose()
        * dataset.std_mag_cell.cell.lattice.basis
        * (magnetic_cell.cell.lattice.basis * dataset.std_linear)
            .try_inverse()
            .unwrap();
    assert_relative_eq!(stretch, stretch.transpose(), epsilon = 1e-12);
    assert!(stretch.symmetric_eigen().eigenvalues.min() > 0.0);
    assert_relative_eq!(
        dataset.std_rotation_matrix.transpose() * dataset.std_rotation_matrix,
        Matrix3::identity(),
        epsilon = 1e-12
    );
    assert_relative_eq!(
        dataset.std_rotation_matrix
            * stretch
            * magnetic_cell.cell.lattice.basis
            * dataset.prim_std_linear,
        dataset.prim_std_mag_cell.cell.lattice.basis,
        epsilon = 1e-10
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
        let dataset = assert_magnetic_dataset_with_default(&magnetic_cell, symprec, action);

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
        let dataset = assert_magnetic_dataset_with_default(&magnetic_cell, symprec, action);

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
        let dataset = assert_magnetic_dataset_with_default(&magnetic_cell, symprec, action);

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
    let action = RotationMagneticMomentAction::Polar;

    let dataset = assert_magnetic_dataset_with_default(&magnetic_cell, symprec, action);

    assert_eq!(dataset.uni_number, 932);
}

#[test]
fn test_with_pyrochlore() {
    let path = Path::new("tests/assets/pyrochlore.json");
    let magnetic_cell: MagneticCell<NonCollinear> =
        serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();

    let symprec = 1e-4;
    let action = RotationMagneticMomentAction::Axial;

    let _dataset = assert_magnetic_dataset_with_default(&magnetic_cell, symprec, action);
}

#[test]
fn test_strained_noncollinear_standardization() {
    let original: MagneticCell<NonCollinear> =
        serde_json::from_str(&fs::read_to_string("tests/assets/pyrochlore.json").unwrap()).unwrap();
    let action = RotationMagneticMomentAction::Axial;
    let reference = MoyoMagneticDataset::with_default(&original, 1e-4, action).unwrap();
    let strain = matrix![
        1.0007, 0.0002, -0.0003;
        0.0002, 0.9995, 0.0001;
        -0.0003, 0.0001, 1.0003;
    ];
    for handedness in [1.0, -1.0] {
        let frame = Rotation3::from_euler_angles(0.37, -0.21, 0.13).into_inner()
            * Matrix3::from_diagonal(&vector![1.0, 1.0, handedness]);
        let input = MagneticCell::new(
            original.cell.lattice.rotate(&(frame * strain)),
            original.cell.positions.clone(),
            original.cell.numbers.clone(),
            original
                .magnetic_moments
                .iter()
                .map(|m| m.act_rotation(&frame, action))
                .collect(),
        );
        let mut outputs = Vec::new();
        for rotate_basis in [false, true] {
            let dataset = MoyoMagneticDataset::new(
                &input,
                0.05,
                AngleTolerance::default(),
                Some(0.02),
                action,
                rotate_basis,
            )
            .unwrap();
            assert_eq!(dataset.uni_number, reference.uni_number);
            for (cell, primitive) in [
                (&dataset.std_mag_cell, false),
                (&dataset.prim_std_mag_cell, true),
            ] {
                let operations =
                    magnetic_operations_from_uni_number(dataset.uni_number, primitive).unwrap();
                let metric = cell.cell.lattice.metric_tensor();
                for operation in operations.iter() {
                    let rotation = operation.operation.rotation.map(f64::from);
                    assert_relative_eq!(
                        rotation.transpose() * metric * rotation,
                        metric,
                        epsilon = 1e-10
                    );
                    let cart_rotation = operation.operation.cartesian_rotation(&cell.cell.lattice);
                    assert_relative_eq!(
                        cart_rotation.transpose() * cart_rotation,
                        Matrix3::identity(),
                        epsilon = 1e-12
                    );
                    for (i, position) in cell.cell.positions.iter().enumerate() {
                        let image = rotation * position + operation.operation.translation;
                        let j = cell
                            .cell
                            .positions
                            .iter()
                            .enumerate()
                            .find_map(|(j, target)| {
                                let delta = image - target;
                                (cell.cell.numbers[i] == cell.cell.numbers[j]
                                    && (delta - delta.map(f64::round)).norm() < 1e-12)
                                    .then_some(j)
                            })
                            .unwrap();
                        let image_moment = cell.magnetic_moments[i].act_magnetic_operation(
                            &cart_rotation,
                            operation.time_reversal,
                            action,
                        );
                        assert_relative_eq!(
                            image_moment.0,
                            cell.magnetic_moments[j].0,
                            epsilon = 1e-12
                        );
                    }
                }
            }
            outputs.push(dataset);
        }
        let [fixed, rotated] = outputs.as_slice() else {
            unreachable!()
        };
        assert_relative_eq!(
            fixed.std_rotation_matrix,
            Matrix3::identity(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            rotated.std_rotation_matrix.determinant(),
            1.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            rotated.std_rotation_matrix * fixed.std_mag_cell.cell.lattice.basis,
            rotated.std_mag_cell.cell.lattice.basis,
            epsilon = 1e-12,
        );
        for (fixed_moment, rotated_moment) in fixed
            .std_mag_cell
            .magnetic_moments
            .iter()
            .zip(&rotated.std_mag_cell.magnetic_moments)
        {
            assert_relative_eq!(
                rotated.std_rotation_matrix * fixed_moment.0,
                rotated_moment.0,
                epsilon = 1e-12
            );
        }
    }
}

#[test]
fn test_strained_left_handed_monoclinic_standardization() {
    let action = RotationMagneticMomentAction::Axial;
    let frame = Rotation3::from_euler_angles(0.37, -0.21, 0.13).into_inner();
    let make_cell = |basis| {
        MagneticCell::new(
            Lattice::from_basis(basis).rotate(&frame),
            vec![vector![0.0, 0.0, 0.0]],
            vec![1],
            vec![NonCollinear(frame * vector![0.0, 1.0, 0.0])],
        )
    };
    let original = make_cell([[4.0, 0.0, 0.0], [0.0, 5.0, 0.0], [-0.7, 0.0, -6.0]]);
    let reference = MoyoMagneticDataset::with_default(&original, 1e-4, action).unwrap();
    let input = make_cell([
        [4.01, 0.002, 0.0],
        [0.003, 5.01, 0.002],
        [-0.7, 0.004, -6.01],
    ]);

    // A twofold axis along b removes the ab and bc metric terms, preserving ac.
    let mut expected_metric = input.cell.lattice.metric_tensor();
    expected_metric[(0, 1)] = 0.0;
    expected_metric[(1, 0)] = 0.0;
    expected_metric[(1, 2)] = 0.0;
    expected_metric[(2, 1)] = 0.0;
    assert!(expected_metric[(0, 2)].abs() > 1.0);

    for rotate_basis in [false, true] {
        let dataset = MoyoMagneticDataset::new(
            &input,
            0.05,
            AngleTolerance::default(),
            Some(0.02),
            action,
            rotate_basis,
        )
        .unwrap();
        assert_eq!(dataset.uni_number, reference.uni_number);
        for (cell, linear) in [
            (&dataset.std_mag_cell, dataset.std_linear),
            (&dataset.prim_std_mag_cell, dataset.prim_std_linear),
        ] {
            assert!(cell.cell.lattice.basis.determinant() < 0.0);
            assert_relative_eq!(
                cell.cell.lattice.metric_tensor(),
                linear.transpose() * expected_metric * linear,
                epsilon = 1e-10
            );
        }
    }
}

#[test]
fn test_with_large_mag_symprec() {
    // https://github.com/spglib/moyo/issues/295
    // With these borderline tolerances the magnetic operation set may not be closed.
    // The requirement is only that `MoyoMagneticDataset::new` does not panic; it may
    // either succeed or return a `MoyoError` such as `TooLargeToleranceError`.
    let magnetic_cell = MagneticCell::new(
        Lattice::from_basis([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        vec![[0.0, 0.0, 0.0].into(), [0.25, 0.25, 0.25].into()],
        vec![1, 1],
        vec![Collinear(0.003), Collinear(0.005)],
    );
    let symprec = 1e-2;
    let angle_tolerance = AngleTolerance::default();
    let mag_symprec = 1e-2;
    let action = RotationMagneticMomentAction::Axial;
    let rotate_basis = false;

    let _ = MoyoMagneticDataset::new(
        &magnetic_cell,
        symprec,
        angle_tolerance,
        Some(mag_symprec),
        action,
        rotate_basis,
    );
}
