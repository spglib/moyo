# Change Log

All notable changes to this project will be documented in this file.
Entries are grouped by workspace crate (`moyo`, `moyopy`, `moyoc`, `moyo-wasm`).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unlisted patch versions (e.g. v0.7.1, v0.7.3, v0.7.6, v0.7.7) contain only dependency or build updates with no user-visible changes.

## Unreleased

### moyo

- Add `hall_symbol()` and `arithmetic_crystal_class()` to `MoyoDataset` (#272)
- Fix panic in `MoyoMagneticDataset::new` when the magnetic operation set is not closed; return `MoyoError::TooLargeToleranceError` instead (#295)

### moyopy

- Add Python binding for `integral_normalizer` (#256)

## v0.7.8 - 2026-03-03

### moyopy

- Add Python 3.14 support (#247)

## v0.7.5 - 2026-02-07

### moyo

- Make `std_cell` idempotent (#235)

## v0.7.4 - 2025-12-26

### moyo

- Fix symmetrization for P-1 (#206)

## v0.7.2 - 2025-11-17

### moyo

- Fix handedness in `std_rotation_matrix` (#194)

## v0.7.0 - 2025-10-25

### moyo

- Change default setting into ITA's standard (#178)
- Add flag to rotate basis in standardization (#180)

### moyopy

- Drop Python 3.9 and support Python 3.13 (#177)
