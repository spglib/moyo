use std::cell::Cell as WarningCount;
use std::sync::Once;

use log::{Level, LevelFilter, Log, Metadata, Record};
use nalgebra::vector;
use rstest::rstest;

use moyo::base::{
    AngleTolerance, Cell, Collinear, Lattice, MagneticCell, Operation, RotationMagneticMomentAction,
};
use moyo::data::{LayerSetting, Setting};
use moyo::{MoyoDataset, MoyoLayerDataset, MoyoMagneticDataset};

struct WarningLogger;

thread_local! {
    static WARNINGS: WarningCount<usize> = const { WarningCount::new(0) };
}

impl Log for WarningLogger {
    fn enabled(&self, metadata: &Metadata<'_>) -> bool {
        metadata.level() == Level::Warn
    }

    fn log(&self, record: &Record<'_>) {
        if self.enabled(record.metadata()) {
            assert!(record.args().to_string().contains("non-integer rotation"));
            WARNINGS.with(|count| count.set(count.get() + 1));
        }
    }

    fn flush(&self) {}
}

fn start_warning_capture() {
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        log::set_logger(&WarningLogger).unwrap();
        log::set_max_level(LevelFilter::Warn);
    });
    WARNINGS.with(|count| count.set(0));
}

fn assert_valid_operations(cell: &Cell, operations: &[Operation]) {
    let metric = cell.lattice.metric_tensor();
    for operation in operations {
        let rotation = operation.rotation.map(f64::from);
        assert!((rotation.determinant().abs() - 1.0).abs() < 1e-8);
        assert!(
            (rotation.transpose() * metric * rotation - metric)
                .abs()
                .max()
                < 1e-8
        );

        let mut visited = vec![false; cell.num_atoms()];
        for (position, number) in cell.positions.iter().zip(&cell.numbers) {
            let transformed = rotation * position + operation.translation;
            let matches: Vec<_> = cell
                .positions
                .iter()
                .zip(&cell.numbers)
                .enumerate()
                .filter(|(_, (candidate, species))| {
                    let diff = transformed - *candidate;
                    *species == number && (diff - diff.map(f64::round)).abs().max() < 1e-8
                })
                .map(|(index, _)| index)
                .collect();
            let [index] = matches.as_slice() else {
                panic!("operation must map each atom to exactly one atom of the same species");
            };
            assert!(!visited[*index]);
            visited[*index] = true;
        }
    }
}

fn rocksalt_supercell(n: usize) -> Cell {
    let base = [
        (vector![0.0, 0.0, 0.0], 11),
        (vector![0.5, 0.5, 0.0], 11),
        (vector![0.5, 0.0, 0.5], 11),
        (vector![0.0, 0.5, 0.5], 11),
        (vector![0.5, 0.5, 0.5], 17),
        (vector![0.0, 0.0, 0.5], 17),
        (vector![0.0, 0.5, 0.0], 17),
        (vector![0.5, 0.0, 0.0], 17),
    ];
    let mut positions = vec![];
    let mut numbers = vec![];
    for copy in 0..n {
        for (position, number) in base {
            positions.push(vector![
                (position[0] + copy as f64) / n as f64,
                position[1],
                position[2]
            ]);
            numbers.push(number);
        }
    }
    Cell::new(
        Lattice::from_basis([
            [5.64 * n as f64, 0.0, 0.0],
            [0.0, 5.64, 0.0],
            [0.0, 0.0, 5.64],
        ]),
        positions,
        numbers,
    )
}

#[rstest]
#[case::conventional(1, 192)]
#[case::double(2, 128)]
#[case::triple(3, 192)]
#[case::quadruple(4, 256)]
#[case::quintuple(5, 320)]
fn test_rocksalt_anisotropic_supercell(#[case] n: usize, #[case] expected_operations: usize) {
    let cell = rocksalt_supercell(n);
    start_warning_capture();
    let dataset =
        MoyoDataset::new(&cell, 1e-3, AngleTolerance::Default, Setting::Spglib, true).unwrap();
    assert_eq!(dataset.number, 225);
    assert_eq!(dataset.operations.len(), expected_operations);
    WARNINGS.with(|count| assert_eq!(count.get(), usize::from(n > 1)));
    assert_valid_operations(&cell, &dataset.operations);
}

#[rstest]
#[case::conventional(1)]
#[case::triple(3)]
fn test_magnetic_rocksalt_supercell(#[case] n: usize) {
    let cell = rocksalt_supercell(n);
    let moments = vec![Collinear(0.0); cell.num_atoms()];
    let magnetic_cell = MagneticCell::new(cell.lattice, cell.positions, cell.numbers, moments);
    start_warning_capture();
    let dataset = MoyoMagneticDataset::new(
        &magnetic_cell,
        1e-3,
        AngleTolerance::Default,
        None,
        RotationMagneticMomentAction::Polar,
        true,
    )
    .unwrap();
    // With zero moments, each compatible spatial operation occurs with and
    // without time reversal.
    assert_eq!(dataset.magnetic_operations.len(), 384);
    WARNINGS.with(|count| assert_eq!(count.get(), usize::from(n > 1)));
    let operations: Vec<_> = dataset
        .magnetic_operations
        .into_iter()
        .map(|operation| operation.operation)
        .collect();
    assert_valid_operations(&magnetic_cell.cell, &operations);
}

#[test]
fn test_skew_simple_cubic_supercell() {
    let cell = Cell::new(
        Lattice::from_basis([[1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 1.0, 2.0]]),
        vec![
            vector![0.0, 0.0, 0.0],
            vector![0.0, 5.0 / 6.0, 0.5],
            vector![0.0, 1.0 / 3.0, 0.0],
            vector![0.0, 1.0 / 6.0, 0.5],
            vector![0.0, 2.0 / 3.0, 0.0],
            vector![0.0, 0.5, 0.5],
        ],
        vec![84; 6],
    );
    start_warning_capture();
    let dataset =
        MoyoDataset::new(&cell, 1e-5, AngleTolerance::Default, Setting::Spglib, true).unwrap();
    assert_eq!(dataset.number, 221);
    assert_eq!(dataset.operations.len(), 24);
    WARNINGS.with(|count| assert_eq!(count.get(), 1));
    assert_valid_operations(&cell, &dataset.operations);
}

#[rstest]
#[case::primitive(1, 16)]
#[case::supercell(3, 24)]
fn test_layer_supercell_warning(#[case] n: usize, #[case] expected_operations: usize) {
    let cell = Cell::new(
        Lattice::from_basis([[n as f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 5.0]]),
        (0..n)
            .map(|i| vector![i as f64 / n as f64, 0.0, 0.0])
            .collect(),
        vec![1; n],
    );
    start_warning_capture();
    let dataset = MoyoLayerDataset::new(
        &cell,
        1e-4,
        AngleTolerance::Default,
        LayerSetting::Standard,
        true,
    )
    .unwrap();
    assert_eq!(dataset.number, 61);
    assert_eq!(dataset.operations.len(), expected_operations);
    WARNINGS.with(|count| assert_eq!(count.get(), usize::from(n > 1)));
    assert_valid_operations(&cell, &dataset.operations);
}
