# Layer-Group Architecture

This document captures the design decisions behind moyo's layer-group support: the public surface, the type-level periodicity contract, the crystallographic conventions, the algorithmic structure, and the static databases. Project-management content (milestones, approval status) lives in the design plan; the user-facing standardization contract lives in [`moyopy/docs/layer_standardization.md`](../moyopy/docs/layer_standardization.md).

Reference paper: Fu et al., "Symmetry classification of 2D materials: layer groups versus space groups", *2D Mater.* **11**, 035009 (2024).

Spglib's implementation is referenced for cross-checking, but moyo deliberately does **not** preserve its layer-group conventions — the goal is a cleaner, stricter, simpler design.

## Goals and non-goals

### Goals

- Detect the layer group (one of 80 types) of a 3D crystal structure that has 2D periodicity (a layer system).
- Produce a `MoyoLayerDataset` analogous to `MoyoDataset`: layer-group number, Hall symbol, layer-group operations in the input cell, primitive/standardized layer cells, Wyckoff positions, site symmetry, transformation matrices.
- Reuse moyo's existing primitives (`Cell`, `Operation`, `PrimitiveCell`, `PeriodicKdTree`, identification machinery) rather than forking the pipeline.

### Non-goals (initial release)

- **Magnetic layer groups (1651 magnetic-LG types)** — deferred to a follow-up. The architecture should leave room for it (mirror the `MoyoMagneticDataset` story) but not implement it now.
- **Rod groups, frieze groups, plane groups** — other subperiodic groups; out of scope.
- **Backwards compatibility with spglib's layer-group conventions** (origin choice, axis swap rules, tolerance policies). Diverging cleanly is encouraged where it simplifies the implementation or hardens correctness.
- **Reading layer-group structures from CIF / files.** moyo does not own input parsing — callers construct a `Cell`.

## Public surface

Top-level entry, mirroring `MoyoDataset::new` / `MoyoDataset::with_default`:

```rust
impl MoyoLayerDataset {
    pub fn new(
        cell: &Cell,
        symprec: f64,
        angle_tolerance: AngleTolerance,
        setting: LayerSetting,
        rotate_basis: bool,
    ) -> Result<Self, MoyoError>;

    pub fn with_default(
        cell: &Cell,
        symprec: f64,
    ) -> Result<Self, MoyoError>;
}
```

There is no `AperiodicAxis` parameter. The convention is fixed (see [Crystallographic conventions](#crystallographic-conventions)): `c` is the aperiodic axis and is perpendicular to `a, b`. Inputs that violate the convention are rejected with a descriptive error rather than silently re-oriented. Callers who need to permute basis vectors must do so themselves before constructing a `Cell`; this keeps the input contract unambiguous and the algorithm code free of axis-permutation bookkeeping.

`LayerSetting` selects which Hall-symbol setting to use when an LG has multiple settings (see [Layer-group identification](#layer-group-identification)).

`MoyoLayerDataset` mirrors `MoyoDataset`. Layer-specific points:

- `number: LayerNumber` (1..=80).
- `hall_number: LayerHallNumber` (own table; see [Hall symbol database](#hall-symbol-database)).
- `hm_symbol: String` (paper Table 5 H-M entry).
- `operations: Operations` (LG operations in the input cell).
- `orbits`, `wyckoffs`, `site_symmetry_symbols`.
- `std_cell`, `std_linear`, `std_origin_shift`, `std_rotation_matrix`.
- `prim_std_cell`, `prim_std_linear`, `prim_std_origin_shift`, `mapping_std_prim`.
- `pearson_symbol`: built from the 2D Bravais type (`mp`, `op`, `oc`, `tp`, `hp`) plus the standardized-cell atom count.
- `symprec`, `angle_tolerance`.

## `Cell` / `Lattice` reuse and the periodicity contract

`Cell` and `Lattice` were introduced for the 3D space-group pipeline with an implicit contract: **all three basis vectors are lattice translations**. Layer groups break that contract — the third basis vector is the aperiodic stacking direction, perpendicular to `(a, b)`. Reusing the existing types directly would let a "layer-shaped" input flow into the bulk pipeline (or vice versa) with the type system unable to catch the mistake.

**Decision: introduce two newtypes, `LayerLattice` and `LayerCell`.**

The perpendicularity invariant `c ⊥ a, b` is purely a property of the lattice basis — it does not depend on atoms. It belongs at the `Lattice` layer. The "no spurious translations along `c`" invariant does depend on atoms and stays at the `Cell` layer. Splitting the two newtypes mirrors the existing `Lattice` / `Cell` split.

```rust
pub struct LayerLattice {
    inner: Lattice,
}

impl LayerLattice {
    /// Validates `c ⊥ a, b` and wraps.
    pub fn new(lattice: Lattice, symprec: f64, angle_tolerance: AngleTolerance)
        -> Result<Self, MoyoError>;
    /// Public read-only view of the basis matrix.
    pub fn basis(&self) -> &Matrix3<f64>;
    pub(crate) fn new_unchecked(lattice: Lattice) -> Self;
    // No `lattice() -> &Lattice` accessor, even pub(crate). In-crate sites
    // that need a `&Lattice` (e.g. to call `search_bravais_group`) construct
    // one explicitly via `Lattice { basis: *ll.basis() }`.
}

pub struct LayerCell {
    lattice: LayerLattice,
    positions: Vec<Position>,
    numbers: Vec<AtomicSpecie>,
}

impl LayerCell {
    /// Validates `cell.lattice` against the perpendicularity contract.
    pub fn new(cell: Cell, symprec: f64, angle_tolerance: AngleTolerance)
        -> Result<Self, MoyoError>;
    pub fn lattice(&self) -> &LayerLattice;
    pub fn positions(&self) -> &[Position];
    pub fn numbers(&self) -> &[AtomicSpecie];
    pub fn num_atoms(&self) -> usize;
    pub(crate) fn new_unchecked(lattice: LayerLattice, positions: Vec<Position>, numbers: Vec<AtomicSpecie>) -> Self;
}
```

- `LayerLattice::new` and `LayerCell::new` validate the perpendicularity contract. Construction is the only public way to get either newtype, so any function taking `&LayerLattice` / `&LayerCell` can rely on the invariant.
- **Conversion is one-way, with no in-crate shortcut.** `Lattice -> LayerLattice` and `Cell -> LayerCell` are allowed (via the validating constructors); the reverse is **not** part of the public surface, **and not even available `pub(crate)`**. Neither newtype hands out a borrow of the wrapped bulk type at any visibility. In-crate sites that need a `Cell` reconstruct one as `Cell::new(Lattice { basis: *lc.lattice().basis() }, lc.positions().to_vec(), lc.numbers().to_vec())`. The reconstruction is a one-line conversion at a handful of call sites and makes the layer-to-bulk crossing visible everywhere it occurs. No `into_inner()`, no `Deref`, no `From<LayerCell> for Cell` / `From<LayerLattice> for Lattice`.
- **`LayerCell::lattice()` returns `&LayerLattice`, not `&Lattice`.** This is safe because `LayerLattice` is itself a layer-only type that bulk helpers cannot consume. Public callers needing the bare 3x3 matrix go through `layer_cell.lattice().basis()`.
- All layer-pipeline entry points take the newtypes: `MoyoLayerDataset::new(&LayerCell, ..)`, `LayerPrimitiveCell::new(&LayerCell, ..)`, `LayerPrimitiveSymmetrySearch::new(&LayerCell, ..)`, `search_layer_bravais_group(&LayerLattice, ..)`.
- Plain getters rather than `Deref`: `Deref<Target = Cell>` / `Deref<Target = Lattice>` would let the newtypes silently impersonate the bulk types and re-introduce the conflation we are trying to prevent. Explicit `.lattice()` / `.positions()` keeps the abstraction visible at call sites.

### Why two newtypes rather than one (`LayerCell` only)

Only `LayerCell` is consumed at public entry points — so a single newtype would suffice if the layer-internal code never threaded a bare `Lattice` between modules. In practice the LG-restricted Bravais filter takes a lattice, not a cell, and other in-crate helpers (Minkowski lift, standardization metric refinement) operate on the basis matrix. Keeping those signatures honest — `&LayerLattice` rather than `&Lattice` — is what `LayerLattice` buys us. The cost is one extra ~40-line file; the gain is that the perpendicularity contract is encoded at the same level it is enforced.

### Why not the alternatives

| Option                                  | Why not                                                                                                              |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| Reuse `Cell` / `Lattice` directly       | Type system does not catch bulk/layer conflation; runtime-only validation.                                           |
| `Periodicity` field on `Cell`/`Lattice` | Invasive across the entire bulk pipeline; overkill while rod / frieze groups are out of scope.                       |
| Generic `Cell<P: Periodicity>`          | Type-inference and generic-propagation pain through `MoyoDataset` and the bindings; not justified by the v1 surface. |

If we later add rod or frieze groups, we revisit — a `Periodicity` field becomes more attractive once there are 3+ subperiodic flavors.

### Implications

- `LayerCell::lattice().basis().column(2)` is the aperiodic vector throughout the layer pipeline. It is *not* used as a candidate translation in the primitive-cell search.
- `LayerLattice::new` / `LayerCell::new` are the single boundaries that convert unvalidated bulk types into layer types known to satisfy the contract. Internal helpers downstream may assume the contract holds.
- 3D-shaped operations on `Cell` / `Lattice` (e.g. `Cell::rotate`) do not have a public path through the newtypes. If the layer pipeline needs them, we add forwarder methods (`LayerCell::rotate(&self) -> LayerCell`) so the result stays inside the newtype rather than escaping as a bare `Cell`.

## Crystallographic conventions

These are the places moyo deviates from spglib intentionally. Each is the simpler, less ambiguous choice.

### The aperiodic axis is always `c`

By contract, the third basis vector `c` is the non-lattice direction. This matches the layer-group convention used throughout the paper and in Kopsky-Litvin (ITE): "we take the lattice vectors to be `a` and `b`" and "the non-lattice vector `c_i` is fixed orthogonal to `a_i` and `b_i`" (paper §2). The convention applies uniformly to all 80 LGs.

Where moyo differs from spglib's tolerant input-handling: moyo **rejects** inputs that violate the convention rather than silently permuting the basis. Detection is cheap, the failure mode is clear, and downstream code can assume the convention without re-checking.

Why this is consistent with all 80 LGs (addressing the monoclinic concern):

- The "aperiodic axis" and the "unique symmetry axis" are *different* axes. Confusing the two is the source of the worry that monoclinic LGs might need `b` as aperiodic (because 3D monoclinic *space groups* take `b` as the unique axis). For layer groups this confusion does not apply.
- Monoclinic / oblique LGs (3-7): the 2-fold (or mirror normal) **is** along `c`, which is also the aperiodic axis. No conflict.
- Monoclinic / rectangular LGs (8-18): the 2-fold (or mirror normal) lies in the `ab` plane, along either `a` or `b` depending on Hall-symbol setting (Table 5 axis codes `:a`, `:b`, and analogous `:b-ac` for orthorhombic derivations). The aperiodic axis is *still* `c`. The in-plane label is a separate decision exposed via `LayerSetting`, analogous to the existing `Setting` enum for space groups.
- Orthorhombic and higher: aperiodic axis is `c`; in-plane axes are labeled per Table 5 settings.

Internal benefit: layer-group point-group matrices have the block form `W = [[A, 0], [0, ±1]]` with `A` integer 2x2 and `W_33 = ±1` (paper eq. 4). This block structure is only well-defined relative to a chosen aperiodic axis. Pinning that axis to the third basis vector is the natural and unique-by-construction choice.

### The aperiodic axis is perpendicular to the in-plane axes

Require `c · a = 0` and `c · b = 0` (within `angle_tolerance`). This is the "orthogonal cell" assumption from the paper (lhs of eq. 5 in the paper's discussion of the standard CCS for layer groups).

Rationale:

- Layer-group point operations must preserve the 2D lattice plane, i.e. they are simultaneously block-orthogonal. If `c` has an in-plane component, the same physical operation no longer has the block form `W_i3 = 0`, and the identification logic must un-shear before doing anything else. Spglib bends over backwards to handle the oblique case; the paper actually fixes `c` perpendicular to `a, b` for the same reason.
- Slab calculations in DFT codes universally use a `c` that is perpendicular to the in-plane axes (or the user can trivially make it so). The constraint is not burdensome and rules out a swathe of bug-prone code paths.

The constraint is checked once at entry. If violated, return an error with a clear message.

### Only `p` and `c` centering exist for layer groups

Layer groups admit only oblique-primitive (`mp`), rectangular-primitive (`op`), rectangular-centered (`oc`), square-primitive (`tp`), and hexagonal-primitive (`hp`) Bravais types. There is no `A`, `B`, `I`, `F`, or `R` in 2D. moyo introduces a separate layer-centering concept rather than reusing the existing 3D centering: the 3D centerings encode `c`-direction lattice translations that don't apply to layer systems, so reusing them would invite errors.

## Algorithmic structure

### In-plane lattice reduction

For internal refinement of the in-plane lattice we use Minkowski reduction (not Delaunay). A dedicated 2D Lagrange-Gauss routine (`minkowski_reduce_2d` in `moyo/src/math/minkowski.rs`) operates on the 2x2 in-plane block and leaves `c` untouched. The flow is: extract `(a, b)` as a 2x2 column-matrix → run `minkowski_reduce_2d` → lift the resulting 2x2 unimodular into a 3x3 with `W_33 = 1`, `W_i3 = W_3i = 0` (paper eq. 4) via `lift_2d_to_3d` so the LG block form is preserved.

Why a dedicated 2D routine rather than calling the existing `minkowski_reduce_greedy` with `U2`: the existing function is *typed* as dimension-generic but its inner CVP update uses hardcoded 3D arithmetic, so it does not actually run for `N = U2`. Writing a small Lagrange-Gauss routine (~30 lines) is simpler than retrofitting the 3D code, and it stays close to the textbook 2D algorithm.

Why Minkowski rather than Delaunay: moyo already uses Minkowski for primitive-cell preparation in `PrimitiveCell::new`, so the entire LG pipeline can stay on a single reduction algorithm. Delaunay is conceptually closer to the paper's appendix B framing, but Minkowski produces the same shortest-basis property in 2D. See the [appendix](#appendix-comparison-of-2d-lattice-reductions) for a side-by-side comparison.

Monoclinic-rectangular note (LG 8-18): the paper's appendix B carves out a special "(b, c) partially Delaunay reduced" rule for this case because it allows inputs where `c` is not perpendicular to the lattice plane. Under our perpendicularity constraint, that concern disappears: combined with the rectangular Bravais constraint (`a ⊥ b`), the in-plane lattice has no shear, so 2D Minkowski reduction is trivial (length sort + sign fix). Which of `a`, `b` is the unique symmetry axis is a labeling decision, not a reduction one — handled by `LayerSetting`.

### Primitive 2D cell search

Analogous to `PrimitiveCell::new`, but candidate translations are filtered to those with the `c`-component (in fractional coords) equal to 0 (within `symprec / |c|`). All other machinery (`PeriodicKdTree`, correspondence solver) carries over.

If a candidate translation has nonzero `c`-component, moyo **rejects the input** rather than try to interpret it — a 2D-periodic structure should not have spurious "stacking" translations along `c`. Concretely: if the user passes an AA-stacked bilayer with `c` being the stacking direction, we want them to either thicken the cell to a true 2D unit or pass `c` as the stacking length; a "translation" of `c/2` along `c` would falsify the layer-group hypothesis.

### Bravais group restricted to LG form

Reuse `search_bravais_group` but post-filter to rotations satisfying `W_i3 = W_3i = 0` (i = 1, 2) and `W_33 = ±1` (paper eq. 4). This excludes cubic point groups and any in-plane / out-of-plane mixing automatically.

Rationale for post-filter rather than re-derivation: the existing search already enumerates all integer rotations preserving the 3D metric tensor; we just keep the layer-compatible ones. Cheap, robust, and shares all the tolerance bookkeeping that has already been debugged.

### Layer-group identification

Mirrors `SpaceGroup::new`:

1. From primitive layer operations, project rotations and run the bulk `PointGroup::new` to get the **geometric crystal class** (only). The bulk routine's `arithmetic_number` and `prim_trans_mat` are *not* reused: the bulk normalization re-orients the layer's `c`-axis (a layer `pm2m` and a layer `pmm2` both project to the same 3D arithmetic class but differ by which in-plane axis carries the 2-fold).
1. For each candidate Hall number returned by `LayerSetting::hall_numbers()`, filter by matching geometric crystal class via the layer arithmetic class (`layer_arithmetic_crystal_class_entry`).
1. For each survivor, sweep a small layer-valid normalizer set (`LAYER_CORRECTION_MATRICES`: identity and `b a -c`, the in-plane `a` ↔ `b` swap with `c` sign-flip) and call `match_origin_shift`.
1. Return the matched LG number, Hall number, and `UnimodularTransformation`.

`LayerGroup::new` lives at `moyo/src/identify/layer_group.rs`.

#### Hall symbol parser

`LayerHallSymbol` (in `moyo/src/data/hall_symbol.rs`) is a thin wrapper that uppercases the lattice prefix (lowercase `p`/`c` is paper convention) and delegates to the bulk `HallSymbol`. `Centering::{P,C}` and `LayerCentering::{P,C}` agree on `linear()` and `lattice_points()`, so traverse / primitive-projection results are layer-correct directly.

#### `LayerSetting`

Variants mirror the bulk `Setting`:

- **`Standard`** — BCS / ITE convention. Sources:
  - **Cell choice 1** for the two monoclinic-oblique LGs with multiple cell choices in ITE — LG 5 (`p11a`) and LG 7 (`p112/a`) — per de la Flor, Souvignier, Madariaga & Aroyo, *Acta Cryst.* **A77**, 559-571 (2021), §2 (i).
  - **Origin choice 2** for centrosymmetric LGs with two ITE origins — LG 52 (`p4/n`), LG 62 (`p4/nbm`), LG 64 (`p4/nmm`) — per the same paper, §2 (ii).
  - **Monoclinic-rectangular `:a`** axis labelling (LG 8-18) per Fu et al. 2024 Table 5 / ITE (which only depicts the `:a` diagram for these LGs). Not specified by the BCS paper.
  - **Orthorhombic `abc`** (LG 19-48) per ITE default labelling.
- **`Spglib`** — smallest `LayerHallNumber` for each LG (spglib's first row per LG). Differs from `Standard` only at LG 52, 62, 64 (origin choice 1 vs 2).
- **`HallNumber(LayerHallNumber)`** — explicit user override.

`Default` for `LayerSetting` is `Standard`.

#### Known v1 gap

The basis-correction sweep covers only the `p`-centered axis-swap cases (monoclinic-rectangular `:a` ↔ `:b` and orthorhombic `abc` ↔ `b-ac`). For `c`-centered LGs (10, 13, 18, 26, 35, 36, 47, 48) the conjugation of `b a -c` by the centering reduces to a group-internal rotation in the primitive basis; alignment for those needs a centering-aware correction matrix. Tracked as future work; existing round-trip tests pin specific Hall numbers to avoid exercising this path.

### Standardization

Mirrors `StandardizedCell::new`. The user-facing output contract lives in [`moyopy/docs/layer_standardization.md`](../moyopy/docs/layer_standardization.md); this section describes how it is produced. Differences from the space-group standardization:

- Use the `LayerCentering` type for centering transforms (no body/face/etc.).
- For triclinic and monoclinic-oblique LGs, run 2D Minkowski reduction on the in-plane `(a, b)` with `c` left invariant. For other LG crystal systems, follow paper Figure 1 / Appendix C metric conditions.
- The c-axis is rotated to Cartesian z when `rotate_basis = true`; in-plane axes are then rotated per the user-doc convention (e.g. `a` along Cartesian x for orthorhombic).
- `|c|` is preserved from the input.

### Wyckoff positions

Source: spglib's `database/layer_Wyckoff.csv` covers all 80 LGs in a parsable CSV format. moyo vendors the CSV and parses it into the same `WyckoffPositionSpace` representation already used for space groups. Upstream traces back to Kopsky-Litvin *International Tables for Crystallography vol. E*.

## Layer-group databases

The LG implementation needs five static databases, mirroring moyo's existing space-group tables in `data/`. Each is one-time data ingest with cross-checking against an independent source. All five live alongside the existing space-group tables (separate files, parallel naming). The implementation prefers extending the existing data layer (e.g. reusing `GeometricCrystalClass`) over duplicating types.

### Geometric / arithmetic crystal classes

- **10 geometric crystal classes** for LGs (paper Table 2): subset of the existing 32 `GeometricCrystalClass` enums, dropping the cubic ones. Reuse the existing enum directly — no new type needed.
- **43 arithmetic-geometric crystal classes** for LGs (paper Table 4): pairs of (geometric class, layer Bravais type), with orientation-distinct variants like `pmm2` vs `pm2m` listed as separate entries — the same convention moyo's existing `ArithmeticCrystalClassEntry` uses for 3D (73 entries, where e.g. `mm2C` (#14) and `2mmC` (#15) are split). Verified against Fu et al. Table 4 (43 rows) and the IUCr Dictionary's definition of arithmetic crystal class. Numbering: 1..=43, ordered as in paper Table 4.

### Lattice systems and Bravais types

- **4 layer lattice systems** (paper Table 3): oblique, rectangular, square, hexagonal. Add a `LayerLatticeSystem` enum (the 3D `LatticeSystem` does not apply — triclinic vs. monoclinic is a *crystal-system* distinction, not a *lattice-system* one for layers).
- **5 layer Bravais types**: `mp`, `op`, `oc`, `tp`, `hp`. Add a `LayerBravaisClass` enum.

### Layer centering

**2 variants**: primitive (`p`) and centered (`c`). `LayerCentering` enum with methods analogous to `Centering::linear()` and `Centering::lattice_points()`. Centering vectors are 2D-only (no `c`-direction lattice translations).

### Hall symbol database

- **116 entries** ingested from spglib's `database/layer_spg.csv`, pinned to commit `12355c77fb7c505a55f52cae36341d73b781a065`. The CSV is row-aligned with paper Fu et al. 2024 Table 5.
- Schema mirrors `HallSymbolEntry`: `hall_number`, `number` (LG 1-80), `arithmetic_number` (1-43), `setting`, `hall_symbol`, `hm_short`, `hm_full`, `centering: LayerCentering`.
- `LayerHallNumber` starts at 1 (spglib's CSV order). The first row for each LG number is the smallest Hall (== `LayerSetting::Spglib` default).
- Lattice symbol convention: lowercase `p`/`c` to distinguish from space-group `P`/`C` (paper convention).
- Ingest: `scripts/layer_hall_symbol.py` fetches `layer_spg.csv` over HTTP from a pinned spglib commit and emits the `LAYER_HALL_SYMBOL_DATABASE` Rust source. The LG → layer-arithmetic-class mapping is hard-coded in the script as 43 contiguous LG ranges, cross-checked against the `// LG <n>` smallest-LG comments in `moyo/src/data/layer_arithmetic_crystal_class.rs`.

### Wyckoff position database

- **628 entries** ingested from spglib's `database/layer_Wyckoff.csv` (same pinned commit). Keyed by `LayerHallNumber` (1-116), so the numbering carries over from the Hall DB directly.
- Schema (`LayerWyckoffPosition`): `hall_number`, `multiplicity`, `letter`, `site_symmetry`, `coordinates` (representative `coord1` only — the orbit is recovered by integer-offset search inside `assign_layer_wyckoff_position`, mirroring the bulk pipeline). The reused `WyckoffPositionSpace` parser consumes `coordinates` unchanged.
- Ingest: `scripts/layer_wyckoff.py` fetches the CSV and emits the `LAYER_WYCKOFF_DATABASE` Rust source.

## Out of scope

- **Detecting layer-group symmetry of a "thick slab" embedded in a 3D periodic cell with vacuum padding.** The caller is responsible for passing a structure that genuinely has 2D translational symmetry. (Explored as a follow-up; the natural API is "given a 3D structure with `c` much larger than vacuum threshold, slice and re-detect".)
- **"Symmorphic representative" / space-group fallback** when the LG cannot be matched. moyo returns an error. Spglib has heuristic fallbacks; we don't.
- **Frieze, rod, plane groups.**

## Risk register

| Risk                                                | Mitigation                                                                                                                                             |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Wyckoff data ingest is tedious and error-prone      | Cross-check against Bilbao web-server output for 80 LGs in a one-shot script; commit a fixture comparing all positions.                                |
| Paper Table 5 Hall symbols may have typos           | Verify each entry by expanding it and checking generators against an independent source (Bilbao). Regenerate the table programmatically once verified. |
| Tolerance handling diverges between LG and SG paths | Reuse `iterative_symmetry_search` rhythm and the same `(symprec, angle_tolerance)` tolerance API; do not invent new parameters.                        |
| Users with non-orthogonal `c` will be unhappy       | Document loudly; provide a helper `Cell::orthogonalize_aperiodic_axis` that warns and returns a corrected cell. Out of scope for v1 if costly.         |

## Validation strategy

The payoff of building this fresh (rather than mirroring spglib quirks) is that moyo can verify against multiple authorities:

1. **Round-trip on database generators.** For each LG Hall number, expand the Hall symbol to operations, build a synthetic `Cell` from the generators acting on a single dummy atom, run identification, assert the identified Hall number matches.
1. **Port spglib's layer-group tests.** Vendor spglib's layer-group test cases as Rust integration tests in moyo. Each ported case becomes an `rstest`-parameterized test asserting (a) identified LG number, (b) Hall number under each `LayerSetting` exposed, and (c) Wyckoff letters per atom.
1. **Cross-check against spglib at runtime.** For inputs that satisfy our stricter preconditions (orthogonal aperiodic axis along `c`), our LG number must equal spglib's. Implement as an optional dev test gated on a feature flag.
1. **Property tests** for 2D Minkowski reduction: any reduced basis is reduced under repeated application; the reduction is invariant under input permutation.

## Appendix: Comparison of 2D lattice reductions

Reference for the [in-plane lattice reduction](#in-plane-lattice-reduction) decision to use 2D Minkowski.

|                                | 2D Minkowski                                                                                                                                         | 2D Niggli                                       | 2D Delaunay                                   |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- | --------------------------------------------- |
| Equivalent classical algorithm | Lagrange-Gauss greedy                                                                                                                                | Lagrange-Gauss + sign-fixing rules              | Selling reduction (superbase)                 |
| Resulting interior angle       | `[60°, 90°]` (with `a·b ≥ 0`) or `[60°, 120°]` (signed)                                                                                              | `[60°, 90°]` after sign fixing                  | `[90°, 120°]`                                 |
| Shortest-basis property        | yes (the two shortest linearly independent vectors)                                                                                                  | yes (same as Minkowski)                         | not strictly, but the superbase contains them |
| Output uniqueness              | up to ±a, ±b, swap                                                                                                                                   | unique (sign rules pin the choice)              | up to ±a, ±b, swap                            |
| Iteration count                | O(log)                                                                                                                                               | O(log)                                          | O(log)                                        |
| In moyo today                  | dedicated 2D `minkowski_reduce_2d`; 3D `minkowski_reduce_greedy` is typed-generic but has hardcoded 3D arithmetic so does not actually work for `U2` | 3D-only (`niggli_reduce`)                       | 3D-only (`delaunay_reduce`)                   |
| Used by Fu et al. (paper)      | --                                                                                                                                                   | --                                              | yes (appendix B)                              |
| Implementation effort for LG   | low — ~30 lines of 2D Lagrange-Gauss                                                                                                                 | medium — would require porting sign rules to 2D | low — but adds a separate algorithm           |

Key observations:

- In 2D, Minkowski and Niggli are essentially the same algorithm (Lagrange-Gauss). Niggli adds sign-fixing conventions on top, which moyo does not need at the reduction step (Hall-symbol settings via `LayerSetting` handle labeling).
- Delaunay produces a "fatter" cell (angles in `[90°, 120°]`) while Minkowski / Niggli produce "tighter" cells (angles in `[60°, 90°]`). For hexagonal lattices both happen to give the canonical 120° cell. For oblique triclinic, they differ cosmetically but encode the same lattice.
- Downstream identification (point group, arithmetic class, Hall match) does not depend on which of the three is used. All three produce equivalent shortest-basis information; the standardization step re-imposes per-LG metric conditions anyway.
- moyo's 3D Minkowski cannot be reused as-is for 2D (hardcoded 3D arithmetic in the inner CVP update). A dedicated 2D Lagrange-Gauss routine is small (~30 lines).

**Decision: stay with Minkowski.** Same algorithm family as the rest of the pipeline (`PrimitiveCell::new`), and the angle-convention difference vs. the paper's Delaunay is cosmetic — not visible at the API level once `LayerSetting` handles axis labeling.
