<div align="center">

# <img src="./logo.svg" alt="moyo" width=450>

[![image](https://img.shields.io/pypi/l/moyopy.svg)](https://pypi.python.org/pypi/moyopy)
[![moyo at crates.io](https://img.shields.io/crates/v/moyo.svg)](https://img.shields.io/crates/v/moyo)
[![image](https://img.shields.io/pypi/v/moyopy.svg)](https://pypi.python.org/pypi/moyopy)

</div>

A fast and robust crystal symmetry finder, written in Rust.

<figure>
    <figcaption><a href="bench/mp/analysis.ipynb">Several times faster symmetry detection than Spglib for Materials Project dataset</a></figcaption>
</figure>

## Interfaces

Each interface has its own README with installation, usage, and development
instructions:

- [Rust (`moyo`)](moyo/README.md): core implementation
  ([crates.io](https://crates.io/crates/moyo),
  [docs](https://docs.rs/moyo/latest/moyo/))
- [Python (`moyopy`)](moyopy/README.md): Python binding
  ([PyPI](https://pypi.org/project/moyopy/),
  [docs](https://spglib.github.io/moyo/python/))
- [C (`moyoc`)](moyoc/README.md): C binding, built and consumed with CMake
  ([docs](https://spglib.github.io/moyo/c/))
- [Fortran](moyoc/README.md#fortran-interface): Fortran interface on top of
  the C binding (optional `moyo::moyof` CMake target)
- [JavaScript (`moyo-wasm`)](moyo-wasm/README.md): JavaScript and WebAssembly
  binding ([npm](https://www.npmjs.com/package/@spglib/moyo-wasm))
- [Web viewer](apps/web/README.md): browse the 230 space groups, 80 layer
  groups, and 1651 magnetic space groups at
  <https://spglib.github.io/moyo/viewer/> (static SPA powered by `moyo-wasm`)

A landing hub linking all interfaces is served at the site root,
<https://spglib.github.io/moyo/> (source in [`apps/landing`](apps/landing)).

The functionality categories below mirror the sections of the [moyopy API
reference](https://spglib.github.io/moyo/python/api.html).

Legend: ✅ supported · 🟡 partial · ❌ not exposed · ➖ not applicable.

| Category             | Functionality        | Rust (`moyo`) | Python (`moyopy`) | C (`moyoc`) | Fortran | JS/WASM (`moyo-wasm`) |
| -------------------- | -------------------- | ------------- | ----------------- | ----------- | ------- | --------------------- |
| Shared               | Core types           | ✅            | ✅                | 🟡          | 🟡      | 🟡                    |
| Space group          | Dataset from cell    | ✅            | ✅                | ✅          | ✅      | ✅                    |
| Space group          | Data access          | ✅            | ✅                | 🟡          | 🟡      | 🟡                    |
| Space group          | Group identification | ✅            | ✅                | 🟡          | 🟡      | ❌                    |
| Layer group          | Dataset from cell    | ✅            | ✅                | ✅          | ✅      | ❌                    |
| Layer group          | Data access          | ✅            | ✅                | 🟡          | 🟡      | 🟡                    |
| Layer group          | Group identification | ✅            | ✅                | ✅          | ✅      | ❌                    |
| Magnetic space group | Dataset from cell    | ✅            | ✅                | ✅          | ✅      | ❌                    |
| Magnetic space group | Data access          | ✅            | ✅                | ✅          | ✅      | 🟡                    |
| Magnetic space group | Group identification | ✅            | ✅                | ✅          | ✅      | ❌                    |

Notes:

- Core types are `Cell`, `Operations`, and `UnimodularTransformation`. C
  exposes `Cell`, magnetic cells, `Operations`, and `MagneticOperations`;
  WASM exposes only `Cell` and (implicitly) `Operations` via the dataset
  output. `UnimodularTransformation` is not bound in either.
- "Dataset from cell" is the main symmetry-analysis entry point
  (`MoyoDataset`, `MoyoLayerDataset`, `MoyoCollinearMagneticDataset`,
  `MoyoNonCollinearMagneticDataset`).
- "Data access" covers classification tables and lookup helpers (e.g.
  `Setting`, `Centering`, `HallSymbolEntry`, `SpaceGroupType`,
  `ArithmeticCrystalClass`, `operations_from_number`, `LayerSetting`,
  `LayerCentering`, `LayerHallSymbolEntry`, `LayerGroupType`,
  `LayerArithmeticCrystalClass`, `operations_from_layer_number`,
  `MagneticSpaceGroupType`, `magnetic_operations_from_uni_number`).
- C's "Data access" is 🟡 because the standalone `Centering`,
  `ArithmeticCrystalClass`, and `WyckoffPosition` classes are not bound
  (their content is reachable as strings through the entry and group-type
  lookups); the settings, Hall symbol entries, group types, and
  `operations_from_*` generators are all exposed.
- "Group identification" recovers a group label from a primitive list of
  symmetry operations (`PointGroup` + `SpaceGroup` + `integral_normalizer`,
  `LayerGroup`, `MagneticSpaceGroup`). C's space-group row is 🟡 because
  `integral_normalizer` is not bound; `PointGroup`, `SpaceGroup`,
  `LayerGroup`, and `MagneticSpaceGroup` are all exposed.
- The Fortran interface wraps the complete moyoc API, so its column matches
  the C column one-to-one.

## How to cite moyo and its interfaces

If you use moyo or its interfaces in your work, please cite the following figshare entry.

```text
@misc{moyo,
  author = {Kohei Shinohara},
  title  = {{moyo: A fast and robust crystal symmetry finder, written in Rust}},
  year   = {2026},
  month  = {1},
  doi    = {10.6084/m9.figshare.31081162.v1},
  url    = {https://figshare.com/articles/software/moyo_A_fast_and_robust_crystal_symmetry_finder_written_in_Rust_/31081162},
  note   = {Source code available at \url{https://github.com/spglib/moyo}}
}
```

This citation may be superseded by a peer-reviewed publication in the future.

## Dev

```shell
cargo install cargo-release cargo-edit cargo-deny cargo-semver-checks
```

### How to release

1. Add user-facing changes under `## Unreleased` in the [changelog](./CHANGELOG.md).
1. Review the release pull request that [tagpr](https://github.com/Songmu/tagpr) maintains. For a minor or major bump, add the `tagpr:minor` or `tagpr:major` label and run the **Prepare release** workflow on `main` to refresh the proposal.
1. Merge the release pull request once checks pass. Tagpr pushes the version tag and the [release workflow](.github/workflows/release.yaml) publishes the packages and GitHub release assets.

#### One-time release setup

1. Install a GitHub App on `spglib/moyo` with read and write access to Contents and Pull requests, and read access to Issues. Save its private key as the `TAGPR_APP_PRIVATE_KEY` Actions secret.
1. Configure [crates.io Trusted Publishing](https://crates.io/docs/trusted-publishing) for `moyo` with GitHub repository `spglib/moyo`, workflow `release.yaml`, and environment `release`.
1. Set the `TAGPR_APP_CLIENT_ID` Actions variable to the App's client ID. Then run **Prepare release** on `main` to open the first release pull request.

Until the App is configured, use `cargo set-version --bump patch` and `cargo semver-checks` to prepare a release pull request. After merging it, run `cargo release --execute` from the updated `main` branch.

### Debugging

```shell
RUST_LOG=debug cargo test -- --nocapture
```

## Acknowledgments

We thank Dr. Yusuke Seto for providing the crystallographic database.
We thank Dr. Juan Rodríguez-Carvajal for providing the magnetic space-group database.
