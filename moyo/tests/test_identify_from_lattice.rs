use nalgebra::{Matrix3, matrix, vector};

use moyo::base::{
    AngleTolerance, Lattice, LayerLattice, Operations, Rotations, UnimodularTransformation,
};
use moyo::data::{
    HallSymbol, LayerSetting, Setting, hall_symbol_entry, magnetic_operations_from_uni_number,
    operations_from_layer_number, operations_from_number,
};
use moyo::identify::{LayerGroup, MagneticSpaceGroup, PointGroup, SpaceGroup};

fn project_rotations(operations: &Operations) -> Rotations {
    operations
        .iter()
        .map(|operation| operation.rotation)
        .collect()
}

fn rebases() -> [UnimodularTransformation; 2] {
    [1, -1].map(|sign| {
        UnimodularTransformation::new(
            matrix![0, sign, 2; 1, 0, 1; 0, 0, -1],
            vector![0.13, 0.27, 0.19],
        )
    })
}

fn assert_operations(actual: &Operations, expected: &Operations, periodic_axes: usize) {
    assert_eq!(actual.len(), expected.len());
    for operation in actual {
        assert!(
            expected.iter().any(|reference| {
                let mut diff = operation.translation - reference.translation;
                for i in 0..periodic_axes {
                    diff[i] -= diff[i].round();
                }
                operation.rotation == reference.rotation && diff.norm() < 1e-8
            }),
            "operation {operation:?} missing from database representative"
        );
    }
}

#[test]
fn space_group_witness_from_lattice() {
    // Screw axes, a cubic group, and a centered monoclinic primitive cell.
    for (number, conventional_basis) in [
        (76, Matrix3::from_diagonal(&vector![4.0, 4.0, 6.0])),
        (
            144,
            matrix![4.0, -2.0, 0.0; 0.0, 12.0_f64.sqrt(), 0.0; 0.0, 0.0, 6.0],
        ),
        (212, Matrix3::from_diagonal(&vector![4.0, 4.0, 4.0])),
        (5, matrix![3.0, 0.0, -1.0; 0.0, 4.0, 0.0; 0.0, 0.0, 5.0]),
        (48, Matrix3::from_diagonal(&vector![3.0, 4.0, 5.0])),
    ] {
        let hall_number = Setting::Standard.hall_number(number).unwrap();
        let symbol = HallSymbol::from_hall_number(hall_number).unwrap();
        let lattice = Lattice {
            basis: conventional_basis
                * symbol
                    .centering
                    .linear()
                    .map(f64::from)
                    .try_inverse()
                    .unwrap(),
        };
        let reference = symbol.primitive_traverse();
        for rebase in rebases() {
            let input = rebase.transform_operations(&reference);
            let input_lattice = rebase.transform_lattice(&lattice);
            for setting in [
                Setting::Standard,
                Setting::Spglib,
                Setting::HallNumber(hall_number),
            ] {
                let group =
                    SpaceGroup::from_lattice(&input_lattice, &input, setting, 1e-8).unwrap();
                assert_eq!(group.number, number);
                assert_eq!(group.transformation.determinant(), rebase.determinant());
                let target =
                    operations_from_number(number, Setting::HallNumber(group.hall_number), true)
                        .unwrap();
                assert_operations(
                    &group.transformation.transform_operations(&input),
                    &target,
                    3,
                );
            }
        }
    }
}

#[test]
fn point_group_witness_from_lattice() {
    // These primitive P settings have the same rotations as their arithmetic
    // point-group representatives (C4, C3, and O, respectively).
    for (number, basis) in [
        (76, Matrix3::from_diagonal(&vector![4.0, 4.0, 6.0])),
        (
            144,
            matrix![4.0, -2.0, 0.0; 0.0, 12.0_f64.sqrt(), 0.0; 0.0, 0.0, 6.0],
        ),
        (212, Matrix3::from_diagonal(&vector![4.0, 4.0, 4.0])),
    ] {
        let reference = operations_from_number(number, Setting::Standard, true).unwrap();
        let reference_rotations = project_rotations(&reference);
        let arithmetic_number = hall_symbol_entry(Setting::Standard.hall_number(number).unwrap())
            .unwrap()
            .arithmetic_number;
        for rebase in rebases() {
            let input = rebase.transform_operations(&reference);
            let group = PointGroup::from_lattice(
                &rebase.transform_lattice(&Lattice { basis }),
                &project_rotations(&input),
            )
            .unwrap();
            assert_eq!(group.arithmetic_number, arithmetic_number);
            let witness = UnimodularTransformation::from_linear(group.prim_trans_mat);
            assert_eq!(witness.determinant(), rebase.determinant());
            let actual = project_rotations(&witness.transform_operations(&input));
            assert_eq!(actual.len(), reference_rotations.len());
            assert!(
                actual
                    .iter()
                    .all(|rotation| reference_rotations.contains(rotation))
            );
        }
    }
}

#[test]
fn magnetic_space_group_witness_from_lattice() {
    // P4_1 reference groups of construct types I, II, III, and IV.
    for uni_number in [667, 668, 669, 670] {
        let reference = magnetic_operations_from_uni_number(uni_number, true).unwrap();
        let lattice = Lattice {
            basis: Matrix3::from_diagonal(&vector![4.0, 4.0, 6.0]),
        };
        for rebase in rebases() {
            let input = rebase.transform_magnetic_operations(&reference);
            let group =
                MagneticSpaceGroup::from_lattice(&rebase.transform_lattice(&lattice), &input, 1e-8)
                    .unwrap();
            assert_eq!(group.uni_number, uni_number);
            assert_eq!(group.transformation.determinant(), rebase.determinant());
            let actual = group.transformation.transform_magnetic_operations(&input);
            for time_reversal in [false, true] {
                let spatial_parts = |operations: &moyo::base::MagneticOperations| {
                    operations
                        .iter()
                        .filter(|op| op.time_reversal == time_reversal)
                        .map(|op| op.operation.clone())
                        .collect()
                };
                assert_operations(&spatial_parts(&actual), &spatial_parts(&reference), 3);
            }
        }
    }
}

#[test]
fn layer_group_witness_from_lattice() {
    let reference = operations_from_layer_number(8, LayerSetting::Standard, true).unwrap();
    let rebase = UnimodularTransformation::new(
        matrix![0, -1, 0; 1, 4, 0; 0, 0, 1],
        vector![0.13, 0.27, 0.19],
    );
    let lattice = LayerLattice::new(
        rebase.transform_lattice(&Lattice {
            basis: Matrix3::from_diagonal(&vector![3.0, 4.0, 5.0]),
        }),
        1e-4,
        AngleTolerance::Default,
    )
    .unwrap();
    let input = rebase.transform_operations(&reference);
    let group = LayerGroup::from_lattice(&lattice, &input, LayerSetting::Standard, 1e-8).unwrap();
    assert_eq!(group.number, 8);
    // Only a and b are periodic: a wrong shift along c must not be hidden by wrapping.
    assert_operations(
        &group.transformation.transform_operations(&input),
        &reference,
        2,
    );
}
