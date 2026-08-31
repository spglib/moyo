use criterion::{Criterion, criterion_group, criterion_main};

use std::fs;
use std::path::Path;

use moyo::MoyoDataset;
use moyo::base::{Cell, Lattice};
use nalgebra::{Matrix3, Vector3};

/// Conventional rocksalt cell repeated `n` times along each axis.
///
/// A perfect supercell is the worst case for the primitive-cell reduction: the
/// pure-translation count grows as `4 * n^3`, and every translation becomes an
/// extra column of the matrix handed to `HNF::new`.
fn rocksalt_supercell(n: usize) -> Cell {
    let base_positions = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.5, 0.5),
        Vector3::new(0.5, 0.0, 0.5),
        Vector3::new(0.5, 0.5, 0.0),
        Vector3::new(0.5, 0.5, 0.5),
        Vector3::new(0.5, 0.0, 0.0),
        Vector3::new(0.0, 0.5, 0.0),
        Vector3::new(0.0, 0.0, 0.5),
    ];
    let base_numbers = [0, 0, 0, 0, 1, 1, 1, 1];

    let lattice = Lattice::new(Matrix3::identity() * (5.0 * n as f64));
    let mut positions = vec![];
    let mut numbers = vec![];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let shift = Vector3::new(i as f64, j as f64, k as f64);
                for (position, number) in base_positions.iter().zip(base_numbers.iter()) {
                    positions.push((position + shift) / n as f64);
                    numbers.push(*number);
                }
            }
        }
    }
    Cell::new(lattice, positions, numbers)
}

pub fn benchmark(c: &mut Criterion) {
    let path = Path::new("tests/assets/mp-1201492.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
    let symprec = 1e-4;
    c.bench_function("dataset_clathrate_Si", |b| {
        b.iter(|| MoyoDataset::with_default(&cell, symprec))
    });

    // Tighter than the clathrate case above: the reported supercell slowdown was
    // measured at 1e-5.
    let supercell_symprec = 1e-5;
    for n in [2, 3, 4] {
        let cell = rocksalt_supercell(n);
        c.bench_function(
            &format!(
                "dataset_rocksalt_supercell_{n}x{n}x{n}_{}_atoms",
                cell.num_atoms()
            ),
            |b| b.iter(|| MoyoDataset::with_default(&cell, supercell_symprec)),
        );
    }
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
