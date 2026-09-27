use std::collections::HashSet;

use nalgebra::{Matrix3, Vector3, matrix, vector};
use rstest::rstest;

use moyo::base::{AngleTolerance, Lattice, Operation, Rotation, UnimodularTransformation};
use moyo::data::{HallSymbol, Setting};
use moyo::identify::Normalizer;

fn operation_keys(
    operations: impl IntoIterator<Item = Operation>,
    free_axes: &[usize],
) -> HashSet<(Rotation, [i64; 3])> {
    operations
        .into_iter()
        .map(|op| {
            let translation = std::array::from_fn(|i| {
                if free_axes.contains(&i) {
                    0
                } else {
                    ((op.translation[i] * 1e8).round() as i64).rem_euclid(100_000_000)
                }
            });
            (op.rotation, translation)
        })
        .collect()
}

fn as_operation(transformation: &UnimodularTransformation) -> Operation {
    Operation::new(*transformation.linear(), *transformation.origin_shift())
}

#[rstest]
#[case::p1(1, vector![3.0, 4.0, 5.0], &[0, 1, 2])]
#[case::p_minus_1(2, vector![3.0, 4.0, 5.0], &[])]
#[case::pmm2(25, vector![3.0, 4.0, 5.0], &[2])]
#[case::pmmm(47, vector![3.0, 4.0, 5.0], &[])]
#[case::p41(76, vector![4.0, 4.0, 6.0], &[2])]
#[case::p432(207, vector![4.0, 4.0, 4.0], &[])]
fn normalizer_is_covariant_under_rebasing(
    #[case] number: i32,
    #[case] lengths: Vector3<f64>,
    #[case] free_axes: &[usize],
    #[values(false, true)] preserve_chirality: bool,
) {
    let lattice = Lattice::new(Matrix3::from_diagonal(&lengths));
    let symbol =
        HallSymbol::from_hall_number(Setting::Standard.hall_number(number).unwrap()).unwrap();
    let operations = symbol.primitive_traverse();
    let generators = symbol.primitive_generators();
    let reference = Normalizer::from_lattice(
        &lattice,
        &operations,
        &generators,
        1e-5,
        AngleTolerance::Default,
        preserve_chirality,
    )
    .unwrap();
    // Compare the whole discrete group, not particular coset representatives
    // or generators. Polar translations are arbitrary, so quotient them out.
    let expected = operation_keys(reference.operations().iter().map(as_operation), free_axes);
    for linear in [
        Matrix3::identity(),
        matrix![-1, 0, 0; 0, 1, 0; 0, 0, 1],
        matrix![0, 1, 0; 1, 0, 0; 0, 0, 1],
        matrix![0, -1, 2; 1, 0, 1; 0, 0, -1],
    ] {
        let rebase = UnimodularTransformation::new(linear, vector![0.13, 0.27, 0.19]);
        let input_lattice = rebase.transform_lattice(&lattice);
        let input_operations = rebase.transform_operations(&operations);
        let input_generators = rebase.transform_operations(&generators);
        let normalizer = Normalizer::from_lattice(
            &input_lattice,
            &input_operations,
            &input_generators,
            1e-5,
            AngleTolerance::Default,
            preserve_chirality,
        )
        .unwrap();
        let all_operations = normalizer.operations();
        let actual = operation_keys(
            all_operations
                .iter()
                .map(|t| rebase.inverse().transform_operation(&as_operation(t))),
            free_axes,
        );
        assert_eq!(actual, expected, "SG {number}, rebase {linear:?}");
        let metric = input_lattice.metric_tensor();
        let input_keys = operation_keys(input_operations.clone(), &[]);
        for transformation in &all_operations {
            let p = transformation.linear_as_f64();
            assert!((p.transpose() * metric * p - metric).norm() < 1e-8);
            assert_eq!(
                operation_keys(transformation.transform_operations(&input_operations), &[]),
                input_keys,
            );
            if preserve_chirality || number == 76 {
                assert_eq!(transformation.determinant(), 1);
            }
        }
        // P4_1 is enantiomorphic, so even its unrestricted normalizer is
        // proper. The other fixtures admit improper physical actions.
        if !preserve_chirality && number != 76 {
            assert!(all_operations.iter().any(|t| t.determinant() == -1));
        }
        assert_eq!(
            normalizer.continuous_translation_directions.len(),
            free_axes.len()
        );
        let directions: Vec<_> = normalizer
            .continuous_translation_directions
            .iter()
            .map(|direction| rebase.linear_as_f64() * direction)
            .collect();
        for direction in &directions {
            for i in 0..3 {
                if !free_axes.contains(&i) {
                    assert!(direction[i].abs() < 1e-8);
                }
            }
        }
        let span = Matrix3::from_fn(|i, j| directions.get(j).map_or(0.0, |d| d[i]));
        assert_eq!(span.svd(false, false).rank(1e-8), free_axes.len());
    }
}
