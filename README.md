# moyo

[![moyo at crates.io](https://img.shields.io/crates/v/moyo.svg)](https://img.shields.io/crates/v/moyo)

## Interfaces

- [Rust](moyo/README.md): core implementation
- [Python](moyopy/README.md)

## Goals

- Find symmetry operations of a given crystal structure, identify its crystallographic group, and symmetrize the given crystal structure
- Well-defined tolerance for finding symmetry operations
- No dependency on existing symmetry-finder packages
- Simplify crystal symmetry algorithm by extensively exploiting the group structures of crystallographic groups

## Non-goals

- Crystallographic groups in other than three dimensions
- Matching two similar crystal structures

## Dev

### How to release

1. Increment the version number in `Cargo.toml`
1. `cargo release --execute`
1. `git push origin <tag-version>`

## Acknowledgments

We thank Dr. Yusuke Seto for providing the crystallographic database.
