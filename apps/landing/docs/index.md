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

    [:octicons-package-16: npm](https://www.npmjs.com/package/@spglib/moyo-wasm)

-   :octicons-eye-16: __Web viewer__

    ---

    Browse the 230 space groups, 80 layer groups, and 1651 magnetic space groups
    in your browser. Powered by `moyo-wasm`; no server, no analytics.

    [:octicons-arrow-right-16: Open the viewer](https://spglib.github.io/moyo/viewer/)

</div>

## Feature support

For the per-interface feature-support matrix and the notes behind it, see the
[repository README](https://github.com/spglib/moyo#interfaces).

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
