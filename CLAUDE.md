# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Moyo is a fast, robust crystal symmetry finder written in Rust. It detects crystallographic symmetries in crystal structures, identifies space groups (230 types) and magnetic space groups (1,651 types), and standardizes/symmetrizes structures.

## Repository Structure

This is a Cargo workspace monorepo with four crates:

- `moyo/` - Core Rust library (the main implementation)
- `moyopy/` - Python bindings via PyO3
- `moyoc/` - C bindings via cbindgen
- `moyo-wasm/` - WebAssembly bindings via wasm-bindgen

## Build and Test Commands

Use `just` to run common tasks (run `just` to list all available commands).

```bash
# Rust
just test           # Run Rust tests
just doc            # Open Rust docs

# Python
just py-build       # Build Python bindings (maturin develop)
just py-install     # Install Python dependencies (uv sync)
just py-test        # Run Python tests

# C
just c-build        # Build C bindings
just c-test         # Run C tests

# Linting
just prek           # Run pre-commit checks (prek)
just prek-all       # Run all checks including manual hooks
```

### Running specific tests

```bash
cargo test test_name              # Run a specific Rust test
RUST_LOG=debug cargo test -- --nocapture  # With debug output
```

## Architecture

Module dependency graph within `moyo/`:

```
math <- base <- data <- identify <- standardize <- lib
        ^---- search <--------------|
```

### Key Modules

- **base**: Core data structures (`Cell`, `Lattice`, `Operation`, `MagneticCell`, `Transformation`)
- **math**: Lattice reduction algorithms (Niggli, Delaunay, Minkowski), integer linear algebra (HNF, SNF)
- **search**: Symmetry operation finding using `PeriodicKdTree`, primitive cell identification
- **identify**: Space group and point group identification from symmetry operations
- **data**: Embedded crystallographic databases (Hall symbols, Wyckoff positions, magnetic space groups)
- **symmetrize**: Structure standardization to conventional/primitive cells

### Main Entry Points

- `MoyoDataset::new()` / `MoyoDataset::with_default()` - Analyze crystal symmetry
- `MoyoMagneticDataset::new()` - Analyze magnetic crystal symmetry

### Key Data Structures

- `Cell`: Crystal structure (lattice basis + fractional positions + atomic numbers)
- `MoyoDataset`: Full symmetry analysis result (space group, operations, Wyckoffs, standardized cells)
- `Operation`: Symmetry operation (3x3 rotation matrix + translation vector)

## Testing

Test fixtures are JSON files in `moyo/tests/assets/`. Integration tests use `rstest` for parameterized testing.

## Release Process

```bash
cargo install cargo-release cargo-edit cargo-deny cargo-semver-checks

# Check for breaking changes
cargo semver-checks

# Bump version
cargo set-version --bump patch

# Release (creates git tag)
cargo release --execute
```
