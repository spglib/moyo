# moyo (Rust)

[![CI](https://github.com/spglib/moyo/actions/workflows/ci-rust.yaml/badge.svg)](https://github.com/spglib/moyo/actions/workflows/ci-rust.yaml)
[![moyo at crates.io](https://img.shields.io/crates/v/moyo.svg)](https://img.shields.io/crates/v/moyo)

The core implementation of moyo in Rust

- Crates.io: https://crates.io/crates/moyo
- Document: https://docs.rs/moyo/latest/moyo/

## Module dependency

```
math <- base <- data <- identify <- standardize <- lib
        ^---- search <--------------|
```

## Goals

- Find symmetry operations of a given crystal structure, identify its crystallographic group, and symmetrize the given crystal structure
- Well-defined tolerance for finding symmetry operations
- No dependency on existing symmetry-finder packages
- Simplify crystal symmetry algorithm by extensively exploiting the group structures of crystallographic groups

## Non-goals

- Crystallographic groups in other than three dimensions
- Matching two similar crystal structures
