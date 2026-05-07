# Plan: Layer Group Support in moyo

Reference paper: Fu et al., "Symmetry classification of 2D materials: layer groups versus space groups", *2D Mater.* **11**, 035009 (2024). PDF at `.claude/resources/Fu_2024_2D_Mater._11_035009.pdf`.

The design and conventions live in two checked-in documents:

- [`docs/layer_architecture.md`](../../docs/layer_architecture.md) — design rationale, type-level periodicity contract, crystallographic conventions, algorithmic structure, database layout, risk register, validation strategy, and the 2D-reduction comparison appendix.
- [`moyopy/docs/layer_standardization.md`](../../moyopy/docs/layer_standardization.md) — user-facing standardized-cell convention (output contract for `MoyoLayerDataset.std_cell` / `prim_std_cell`).

This file is the **project-management view**: section approvals and milestone status. Anything that is a design decision should land in the architecture doc; this file should not duplicate it.

## Approval status

Tick when the section is approved. Author edits a `- [ ]` to `- [x]` once agreement is reached.

- [x] §1. Goals & non-goals
- [x] §2. User-facing API
  - [x] §2.1 Entry points and `MoyoLayerDataset`
  - [x] §2.2 Standardized cell convention (now in [`moyopy/docs/layer_standardization.md`](../../moyopy/docs/layer_standardization.md))
  - [x] §2.3 `Cell` / `Lattice` reuse and the periodicity contract
- [x] §3. Crystallographic conventions
  - [x] §3.1 Aperiodic axis is `c`
  - [x] §3.2 Aperiodic axis is perpendicular
  - [x] §3.3 Only `p`/`c` centering
- [x] §4. Algorithmic changes
  - [x] §4.1 In-plane lattice reduction
  - [x] §4.2 Primitive 2D cell search
  - [x] §4.3 Bravais group restricted to LG form
  - [x] §4.4 Layer-group identification
  - [x] §4.5 Standardization
  - [x] §4.6 Wyckoff positions
- [x] §5. Layer-group databases
  - [x] §5.1 Geometric / arithmetic crystal classes
  - [x] §5.2 Lattice systems and Bravais types
  - [x] §5.3 Layer centering
  - [x] §5.4 Layer Hall symbol database
  - [x] §5.5 Layer Wyckoff position database
  - [x] §5.6 Where they live
- [x] §6. Validation strategy
- [x] §7. Python / C / WASM bindings
- [x] §8. Milestones
- [x] §9. Risk register
- [x] §10. Out of scope
- Appendix A. Comparison of 2D lattice reductions (reference, no approval gate)

Section numbering above refers to the corresponding sections in [`docs/layer_architecture.md`](../../docs/layer_architecture.md). Bindings (§7) follow the M5 milestone.

## Milestones

Each milestone produces a passing test suite and is shippable independently.

### M1: Foundations

**Shipped.**

- Layer-centering, layer-lattice-system, and layer arithmetic crystal class enums and tables.
- 2D Minkowski (Lagrange-Gauss) reduction with property tests.
- Aperiodic-axis validator (orthogonality and "c is non-lattice" checks).

### M2: Symmetry search

**Shipped** (PR #304, branch `feat/layer-group-search`).

- 2D primitive cell search.
- LG-restricted Bravais group filter.
- Primitive-layer symmetry search.
- Round-trip on simple LG examples (`p1`, `p2`, `p4mm`) built by hand.

### M3: Hall symbol DB + identification

**Shipped** (PR #305, branch `feat/layer-group-hall-db`).

- [x] 116-entry `LAYER_HALL_SYMBOL_DATABASE` ingested from spglib's `database/layer_spg.csv` (pinned commit `12355c77`); script at `scripts/layer_hall_symbol.py`.
- [x] `LayerHallSymbol` parser (lowercase `p`/`c`); thin wrapper over the bulk `HallSymbol` via lattice-prefix uppercasing.
- [x] `LayerSetting` enum (`HallNumber(n)` / `Spglib` / `Standard`).
- [x] `LayerGroup` identification: `PointGroup`-based geometric-class filter, `LayerSetting::hall_numbers()` iteration, in-plane normalizer sweep (`LAYER_CORRECTION_MATRICES`), `match_origin_shift`.
- [x] Round-trip on all 116 Hall numbers (`test_round_trip_all_layer_hall_numbers`, pinned per-Hall) and on every Standard default (`test_round_trip_default_hall_numbers_via_standard_setting`).
- [x] Non-canonical-basis alignment for `p`-centered LGs verified (`test_normalizer_aligns_non_canonical_basis`, 12 cases).

**Deferred to M5**: centering-aware correction matrices for `c`-centered LGs (10, 13, 18, 26, 35, 36, 47, 48) so that non-canonical `:b` / `b-ac` inputs in centered cells also round-trip under `LayerSetting::Standard`.

### M4: Standardization + Wyckoff

**Shipped** (branch `feat/layer-group-m4-standardize`).

- [x] 628-entry `LAYER_WYCKOFF_DATABASE` ingested from spglib's `database/layer_Wyckoff.csv` (pinned commit `12355c77`); script at `scripts/layer_wyckoff.py`. Keyed by `LayerHallNumber` (1-116).
- [x] `LayerWyckoffPosition` parser reuses `WyckoffPositionSpace`; the representative coordinate (`coord1`) is stored per entry, and the orbit is recovered by integer-offset search inside `assign_layer_wyckoff_position` (mirrors the bulk pipeline).
- [x] `StandardizedLayerCell` (`moyo/src/symmetrize/layer_standardize.rs`): produces the LG-canonical conventional layer cell with `|c|` preserved, `c` along z, `a` along x (when `rotate_basis = true`); applies `LayerCentering::linear()` for `c`-centered LGs; symmetrizes the in-plane metric via Cholesky on the layer-block-averaged metric tensor and pins `g_13 = g_23 = 0` exactly to keep the aperiodic axis decoupled.
- [x] Site-symmetry symbol + Wyckoff letter assignment per atom in the conventional cell.
- [x] Round-trip tests on synthetic LG 1, 3, 55 cells; in-plane skew removal test; multi-atom orbit collapse test.

### M5: Public API + bindings

- `MoyoLayerDataset` and `detect_layer_group`.
- moyopy and moyoc bindings.
- README / docs update.

### M6 (deferred): Magnetic layer groups

Out of scope for v1.
