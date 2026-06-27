---
hide:
  - navigation
  - toc
---

# moyo

A fast and robust crystal symmetry finder, written in Rust. It detects
crystallographic symmetries, identifies the 230 space groups and 1651 magnetic
space groups, and standardizes and symmetrizes structures.

!!! info "The Web Viewer has moved"

    The interactive Web Viewer now lives at
    [**spglib.github.io/moyo/viewer/**](https://spglib.github.io/moyo/viewer/).
    Please update your bookmarks; the old root URL now serves this landing page.

## Interfaces

<div class="grid cards" markdown>

-   :simple-rust: __Rust (`moyo`)__

    ---

    The core implementation.

    [:octicons-package-16: crates.io](https://crates.io/crates/moyo) &middot;
    [:octicons-book-16: docs.rs](https://docs.rs/moyo)

-   :simple-python: __Python (`moyopy`)__

    ---

    Python binding via PyO3.

    [:octicons-package-16: PyPI](https://pypi.org/project/moyopy/) &middot;
    [:octicons-book-16: Docs](https://spglib.github.io/moyo/python/)

-   :simple-c: __C (`moyoc`)__

    ---

    C binding, built and consumed with CMake.

    [:octicons-book-16: Docs](https://spglib.github.io/moyo/c/)

-   :simple-fortran: __Fortran__

    ---

    Fortran interface on top of the C binding.

    [:octicons-book-16: Docs](https://spglib.github.io/moyo/c/fortran-api/)

-   :simple-javascript: __JS / WASM (`moyo-wasm`)__

    ---

    JavaScript and WebAssembly binding.

    [:octicons-package-16: npm](https://www.npmjs.com/package/@spglib/moyo-wasm) &middot;
    [:octicons-globe-16: Web viewer](https://spglib.github.io/moyo/viewer/)

-   :octicons-eye-16: __Web viewer__

    ---

    Browse the 230 space groups, 80 layer groups, and 1651 magnetic space groups
    in your browser. Powered by `moyo-wasm`; no server, no analytics.

    [:octicons-arrow-right-16: Open the viewer](https://spglib.github.io/moyo/viewer/)

</div>

## Feature support

The functionality categories below mirror the sections of the
[moyopy API reference](https://spglib.github.io/moyo/python/api.html).

Legend: :white_check_mark: supported &middot; :yellow_circle: partial &middot;
:x: not exposed &middot; :heavy_minus_sign: not applicable.

| Category             | Functionality        | Rust (`moyo`) | Python (`moyopy`) | C (`moyoc`) | JS/WASM (`moyo-wasm`) |
| -------------------- | -------------------- | ------------- | ----------------- | ----------- | --------------------- |
| Shared               | Core types           | :white_check_mark: | :white_check_mark: | :yellow_circle: | :yellow_circle: |
| Space group          | Dataset from cell    | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| Space group          | Data access          | :white_check_mark: | :white_check_mark: | :yellow_circle: | :yellow_circle: |
| Space group          | Group identification | :white_check_mark: | :white_check_mark: | :yellow_circle: | :x: |
| Layer group          | Dataset from cell    | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| Layer group          | Data access          | :white_check_mark: | :white_check_mark: | :yellow_circle: | :yellow_circle: |
| Layer group          | Group identification | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| Magnetic space group | Dataset from cell    | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |
| Magnetic space group | Data access          | :white_check_mark: | :white_check_mark: | :white_check_mark: | :yellow_circle: |
| Magnetic space group | Group identification | :white_check_mark: | :white_check_mark: | :white_check_mark: | :x: |

See the [repository README](https://github.com/spglib/moyo#interfaces) for the
detailed notes behind this matrix.

## How to cite

If you use moyo or its interfaces in your work, please cite the following figshare
entry:

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
