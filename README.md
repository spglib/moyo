# moyo

[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/spglib/moyo/main.svg)](https://results.pre-commit.ci/latest/github/spglib/moyo/main)
[![image](https://img.shields.io/pypi/l/moyopy.svg)](https://pypi.python.org/pypi/moyopy)
[![moyo at crates.io](https://img.shields.io/crates/v/moyo.svg)](https://img.shields.io/crates/v/moyo)
[![image](https://img.shields.io/pypi/v/moyopy.svg)](https://pypi.python.org/pypi/moyopy)

A fast and robust crystal symmetry finder, written in Rust.

<figure>
    <div style="text-align: center">
        <img src="bench/mp/mp.png" width=50%>
    </div>
    <figcaption><a href="bench/mp/analysis.ipynb">Several times faster symmetry detection than Spglib for Materials Project dataset</a></figcaption>
</figure>

- Rust support available via [crates.io](https://crates.io/crates/moyo): [Docs](https://docs.rs/moyo/latest/moyo/)
- Python support available via [PyPI](https://pypi.org/project/moyopy/): [Docs](https://spglib.github.io/moyo/python/)

## Interfaces

- [Rust](moyo/README.md): core implementation
- [Python](moyopy/README.md): Python binding
- [C](moyoc/README.md): C binding

## Dev

```shell
cargo install cargo-release cargo-edit cargo-deny cargo-semver-checks
```

### How to release

1. `cargo semver-checks` to lint a new release
1. `cargo set-version --bump patch` for patch version increment
1. Write change log and git-commit
1. `cargo release --execute` (If you already release the package in crates.io, run `cargo release --execute --no-publish`)

### Debugging

```shell
RUST_LOG=debug cargo test -- --nocapture
```

## Acknowledgments

We thank Dr. Yusuke Seto for providing the crystallographic database.
We thank Juan Rodríguez-Carvajal for providing the magnetic space-group database.
