# Change Log

All notable changes to this project will be documented in this file.
Entries are grouped by workspace crate (`moyo`, `moyopy`, `moyoc`, `moyo-wasm`).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Unlisted patch versions (e.g. v0.7.1, v0.7.3, v0.7.6, v0.7.7) contain only dependency or build updates with no user-visible changes.

## Unreleased

### moyo

- Add `iter_wyckoff_positions_from_hall_number` and make `WyckoffPosition` public in `moyo::data` to allow looking up all Wyckoff positions of a space group by Hall number without analyzing a structure

### moyo-wasm

- Add `wyckoff_positions` data-access wrapper (and `MoyoWyckoffPosition` type) returning the Wyckoff positions of a space group for the given Hall number

## v0.10.0 - 2026-05-24

### moyo

- Add Euclidean normalizer of space groups: `Normalizer` type and `MoyoDataset::euclidean_normalizer` (#326)
- Expose `operations_from_number`, `operations_from_layer_number`, and `magnetic_operations_from_uni_number` from `moyo::data` so all bindings can reuse them
- Wrap layer-group symmetry search in an iterative tolerance-retry loop (`iterative_layer_symmetry_search`), recovering from `TooLargeToleranceError` raised inside `LayerPrimitiveCell::new` / `LayerPrimitiveSymmetrySearch::new` instead of propagating immediately (#337)
- Fix `WyckoffPositionAssignmentError` for oblique layer groups (LG 5, 7) in near-square cells: the in-plane Minkowski reduction during standardization could rotate a glide off the canonical direction. The reduction is now applied only when it normalizes the canonical operation set, and its origin shift is composed correctly

### moyo-wasm

- Add data-access wrappers for space groups: `operations_from_number`, `hall_symbol_entry`, `space_group_type`, `arithmetic_crystal_class`
- Add data-access wrappers for layer groups: `operations_from_layer_number`, `layer_hall_symbol_entry`, `layer_group_type`, `layer_arithmetic_crystal_class`
- Add data-access wrappers for magnetic space groups: `magnetic_operations_from_uni_number`, `magnetic_hall_symbol_entry`, `magnetic_space_group_type`

## v0.9.0 - 2026-05-09

### moyo

- Add layer-group identification pipeline: `LayerGroup`, `LayerHallSymbol`, `LayerSetting`, and `MoyoLayerDataset` public API (#303, #304, #305, #306, #308)
- Add `LayerGroup::from_lattice` for in-plane lattice reduction (#318)

### moyopy

- Mark Development Status as Beta (#302)
- Add `MoyoLayerDataset` and `LayerSetting` bindings (#310)
- Add `LayerGroupType` bindings and `operations_from_layer_number` (#317)
- Add `LayerGroup` binding for layer-group identification (#318)

## v0.8.0 - 2026-05-01

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
