use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use nalgebra::{matrix, vector};

use moyo::base::{Cell, Lattice, Position};
use moyo::search::{solve_correspondence, solve_correspondence_naive, PeriodicKdTree};

/// O(num_atoms^3)
fn naive(reduced_cell: &Cell) {
    let num_atoms = reduced_cell.num_atoms();
    let symprec = 1e-5;
    for j in 0..num_atoms {
        let translation = reduced_cell.positions[j] - reduced_cell.positions[0];
        let new_positions: Vec<Position> = reduced_cell
            .positions
            .iter()
            .map(|pos| pos + translation)
            .collect();

        solve_correspondence_naive(reduced_cell, &new_positions, symprec);
    }
}

/// O(num_atoms^2 * log(num_atoms))
fn kdtree(reduced_cell: &Cell) {
    let num_atoms = reduced_cell.num_atoms();
    let symprec = 1e-5;
    let pkdtree = PeriodicKdTree::new(reduced_cell, symprec);
    for j in 0..num_atoms {
        let translation = reduced_cell.positions[j] - reduced_cell.positions[0];
        let new_positions: Vec<Position> = reduced_cell
            .positions
            .iter()
            .map(|pos| pos + translation)
            .collect();

        solve_correspondence(&pkdtree, reduced_cell, &new_positions);
    }
}

fn cell_for_benchmark(n: usize) -> Cell {
    let mut positions = vec![];
    let mut numbers = vec![];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                positions.push(vector![
                    i as f64 / n as f64,
                    j as f64 / n as f64,
                    k as f64 / n as f64
                ]);
                numbers.push(0);
            }
        }
    }

    Cell::new(
        Lattice::new(matrix![
            n as f64, 0.0, 0.0;
            0.0, n as f64, 0.0;
            0.0, 0.0, n as f64;
        ]),
        positions,
        numbers,
    )
}

pub fn benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("translation search");
    for n in 1..=8 {
        let cell = cell_for_benchmark(n);
        group.throughput(Throughput::Elements(cell.num_atoms() as u64));
        group.bench_with_input(BenchmarkId::new("naive", n), &cell, |b, cell| {
            b.iter(|| naive(&cell));
        });
        group.bench_with_input(BenchmarkId::new("kdtree", n), &cell, |b, cell| {
            b.iter(|| kdtree(&cell));
        });
    }
    group.finish();
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
