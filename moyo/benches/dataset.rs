use criterion::{Criterion, criterion_group, criterion_main};

use serde_json;
use std::fs;
use std::path::Path;

use moyo::MoyoDataset;
use moyo::base::{AngleTolerance, Cell};
use moyo::data::Setting;

pub fn benchmark(c: &mut Criterion) {
    let path = Path::new("tests/assets/mp-1201492.json");
    let cell: Cell = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
    let symprec = 1e-4;
    c.bench_function("dataset_clathrate_Si", |b| {
        b.iter(|| MoyoDataset::with_default(&cell, symprec))
    });
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
