# Moyo

[![image](https://img.shields.io/pypi/l/moyopy.svg)](https://pypi.python.org/pypi/moyopy)
[![Moyo at crates.io](https://img.shields.io/crates/v/moyo.svg)](https://img.shields.io/crates/v/moyo)
[![image](https://img.shields.io/pypi/v/moyopy.svg)](https://pypi.python.org/pypi/moyopy)

A fast and robust crystal symmetry finder, written in Rust.

<figure>
    <div style="text-align: center">
        <img src="bench/mp/mp.png" width=50%>
    </div>
    <figcaption><a href="bench/mp/analysis.ipynb">Several times faster symmetry detection than Spglib for Materials Project dataset</a></figcaption>
</figure>

- Support most of the symmetry search functionality in Spglib
  - Primitive cell search
  - Symmetry operation search
  - Space-group type identification
  - Wyckoff position assignment
  - Crystal structure symmetrization
- Rust support available via [crates.io](https://crates.io/crates/moyo)
- Python support available via [PyPI](https://pypi.org/project/moyopy/)

## Interfaces

- [Rust](moyo/README.md): core implementation
- [Python](moyopy/README.md)

## Dev

### How to release

1. `cargo set-version --bump patch` for patch version increment
1. Write change log and git-commit
1. `cargo release --execute`

### Debugging

```shell
RUST_LOG=debug cargo test -- --nocapture
```

## Acknowledgments

We thank Dr. Yusuke Seto for providing the crystallographic database.
