---
name: port-rust-docs-to-pyi
description: Port docstrings between PyO3 bindings (`moyopy/src/**/*.rs`) and the Python type stubs (`moyopy/python/moyopy/*.pyi`). Use when adding/editing a class, getter, method, or `#[pyfunction]` so that `__doc__` (Rust `///`), the `.pyi` stub, and the Sphinx API reference all stay in sync.
---

# Port Rust docs to Python stubs

The moyopy bindings have **two** documentation surfaces that must stay synchronized by hand:

1. **Rust `///` doc comments** in `moyopy/src/**/*.rs` -- PyO3 forwards these to Python `__doc__` (visible from `help()`, IDE hover, Sphinx autoapi).
1. **Hand-written `.pyi` stubs** in `moyopy/python/moyopy/*.pyi` -- read by mypy/Pyright/Pylance for type-checker tooltips.

There is **no auto-generation step**. `pyo3-stub-gen` was evaluated and rejected (it forces `module = "<pkg>._<sub>"`, which leaks `_moyopy` into `__module__`). When you change docs in one place, you must mirror them in the other.

## When to use

- Adding a new `#[pyclass]`, `#[getter]`, `#[pymethods]`, or `#[pyfunction]` to `moyopy/src/`.
- Editing an existing `///` doc comment on a binding.
- Editing an existing docstring in a `.pyi` stub.
- Reviewing a PR that adds Python-facing API and the docs feel uneven.

## File mapping

| Rust binding                                   | Python stub                | Public name(s)                                                          |
| ---------------------------------------------- | -------------------------- | ----------------------------------------------------------------------- |
| `moyopy/src/base/cell.rs`                      | `_base.pyi`                | `Cell`                                                                  |
| `moyopy/src/base/magnetic_cell.rs`             | `_base.pyi`                | `CollinearMagneticCell`, `NonCollinearMagneticCell`                     |
| `moyopy/src/base/operation.rs`                 | `_base.pyi`                | `Operations`, `MagneticOperations`, `UnimodularTransformation`          |
| `moyopy/src/base/error.rs`                     | (not exposed in `.pyi`)    | `MoyoError`                                                             |
| `moyopy/src/data/setting.rs`                   | `_data.pyi`                | `Setting`                                                               |
| `moyopy/src/data/centering.rs`                 | `_data.pyi`                | `Centering`                                                             |
| `moyopy/src/data/hall_symbol.rs`               | `_data.pyi`                | `HallSymbolEntry`                                                       |
| `moyopy/src/data/space_group_type.rs`          | `_data.pyi`                | `SpaceGroupType`                                                        |
| `moyopy/src/data/magnetic_space_group_type.rs` | `_data.pyi`                | `MagneticSpaceGroupType`                                                |
| `moyopy/src/data/arithmetic_crystal_class.rs`  | `_data.pyi`                | `ArithmeticCrystalClass`                                                |
| `moyopy/src/data.rs`                           | `_data.pyi`                | `operations_from_number`, `magnetic_operations_from_uni_number`         |
| `moyopy/src/dataset/space_group.rs`            | `_dataset.pyi`             | `MoyoDataset`                                                           |
| `moyopy/src/dataset/magnetic_space_group.rs`   | `_dataset.pyi`             | `MoyoCollinearMagneticDataset`, `MoyoNonCollinearMagneticDataset`       |
| `moyopy/src/identify.rs`                       | `_identify.pyi`            | `PointGroup`, `SpaceGroup`, `MagneticSpaceGroup`, `integral_normalizer` |
| `moyopy/src/lib.rs`                            | `_moyopy.pyi` (re-exports) | (top-level package)                                                     |

`moyopy/python/moyopy/__init__.py` does `from ._moyopy import *`, and `_moyopy.pyi` re-imports symbols from the per-domain `.pyi` files for type-checker visibility.

## Workflow

### Adding or editing a binding item

1. **Write `///` on the Rust side** above the `#[pyclass]` / `#[pymethods]` block, individual `#[getter]`, method, or `#[pyfunction]`. Place it **after** `#[derive(...)]` but **before** `#[pyclass(...)]`/`#[pymethods]`/`#[getter]`. PyO3 forwards it to `__doc__`.
1. **Mirror it in the matching `.pyi` file** as a Python triple-quoted docstring on the same construct.
1. **Verify** with the smoke test in "Verification" below.

### Editing only the .pyi (rarely correct)

Almost always wrong -- if the docstring isn't in the Rust `///`, it doesn't show in `help()`, IDE hover, or Sphinx. Move the prose into the Rust source first, then mirror back to `.pyi`.

### Editing only the Rust /// (also rarely correct)

Type-checker tooltips will go stale. After editing the Rust side, mirror to the matching `.pyi`.

## Conventions

- **Style**: NumPy-style. Class-level docstring first, then per-method/getter/property. Constructors get a `Parameters` section. Avoid Google-style `Args:` -- the project uses NumPy throughout.
- **Inline code**: double backticks (RST/Sphinx style): ` `cell.basis` `, ` `symprec` `. Markdown single backticks render wrong in Sphinx.
- **Code blocks in Rust `///`**: never use indented blocks (rustdoc treats them as Rust and runs them as doctests). Use fenced ```` ```text ```` blocks:
  ````rust
  /// Same transformation convention as the standardized cell:
  ///
  /// ```text
  /// std_cell.basis.T = std_rotation_matrix @ cell.basis.T @ std_linear
  /// ```
  ````
- **Code blocks in `.pyi`**: RST `::` followed by indent, OR fenced -- both render in Sphinx. Match the existing surrounding style in the file you're editing.
- **Cross-references in `.pyi`**: use `:attr:`, `:class:`, `:func:` (Sphinx). In Rust `///`, plain prose is fine -- the Rust doc tree is not used by the Python project.
- **Module attribute** (CRITICAL): keep `#[pyo3(module = "moyopy")]` on every `#[pyclass]` and `#[pyfunction]`. Do **not** change to `"moyopy._moyopy"` -- it leaks the hidden submodule into `__module__`, `repr()`, and pickle. The reviewer rejected this.
- **No `pyo3-stub-gen` re-introduction**: this was tried and reverted (PR #294). It hard-rejects items declared above `module-name` in `pyproject.toml`, with no escape hatch as of v0.22. If you want auto-generation, do it via a different mechanism (mypy stubgen, custom introspection script).

## Verification

After porting, run all of:

```
# rebuild bindings
PYO3_PYTHON=/path/to/moyopy/.venv/bin/python uv run --directory moyopy maturin develop --release --manifest-path Cargo.toml

# pytest
moyopy/.venv/bin/python -m pytest -q moyopy/python/tests

# spot-check __doc__ at runtime
moyopy/.venv/bin/python -c 'import moyopy; print(moyopy.MoyoDataset.std_cell.__doc__)'

# pre-commit
prek run --all-files

# cargo doc-tests pick up bad indented blocks
PYO3_PYTHON=/path/to/.venv/bin/python PYTHONHOME=$(/path/to/.venv/bin/python -c 'import sys; print(sys.base_prefix)') \
  cargo test --manifest-path moyopy/Cargo.toml
```

If `cargo test` fails on a doc-test like `expected one of '!' or '::', found '.'`, you have an indented code block in a `///` -- fence it as ```` ```text ````.

## Pitfalls

- **Indented code in `///`**: rustdoc treats it as Rust and tries to run it as a doctest. Always use ```` ```text ```` fences for non-Rust example blocks.
- **Forgetting to mirror**: the Rust side and the `.pyi` drift silently if you only edit one. Set yourself a checklist: every binding edit touches at least two files (the `.rs` and the matching `.pyi`).
- **Changing `module = "moyopy"`**: do not. See "Conventions" above.
- **Re-adding `pyo3-stub-gen`**: do not (without addressing the parent-module constraint upstream first).
- **`MoyoError` is not in `.pyi`**: it's `#[pyclass]` but not registered in `#[pymodule]` / not in `__all__`. Don't add it to `_base.pyi` unless you also expose it.
- **Re-exports in `_moyopy.pyi`**: when adding a new class, also add it to `__init__.py`'s `__all__` and to `_moyopy.pyi`'s import list and `__all__`.

## Why two surfaces?

PyO3's `///` -> `__doc__` only covers runtime `help()` / hover. Mypy / Pyright read `.pyi` instead and ignore the `__doc__` of the compiled extension. Until pyo3-stub-gen (or an alternative) fits the project's constraints, both surfaces have to be maintained by hand.
