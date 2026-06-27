# API Reference

The public API mirrors the Rust crate layout: `base` (crystal structures and
operation containers), `dataset` (symmetry analysis results), `data`
(crystallographic classification tables), and `identify` (group identification
primitives). The `moyopy.interface` adapters provide conversion helpers to/from
pymatgen and ASE.

This page documents the core types and the space-group API. Layer-group and
magnetic-space-group APIs live on additional pages:

- [Layer Group](api_layer_group.md)
- [Magnetic Space Group](api_magnetic_space_group.md)

## Core types

Lattice + sites, the non-magnetic operation container, and the unimodular
transformation type used across all group families.

::: moyopy.Cell

::: moyopy.Operations

::: moyopy.UnimodularTransformation

## Symmetry datasets

Run a symmetry analysis on a [`moyopy.Cell`][moyopy.Cell] and inspect the result.

::: moyopy.MoyoDataset

## Crystallographic data

Hall symbols, settings, centering, classification tables, and helpers to fetch
operations by ITA number.

::: moyopy.Setting

::: moyopy.Centering

::: moyopy.HallSymbolEntry

::: moyopy.SpaceGroupType

::: moyopy.ArithmeticCrystalClass

::: moyopy.operations_from_number

## Group identification

Identify point groups and space groups from a primitive list of symmetry
operations.

::: moyopy.PointGroup

::: moyopy.SpaceGroup

::: moyopy.integral_normalizer

## Adapters

Convert between [`moyopy.Cell`][moyopy.Cell] and pymatgen `Structure` / ASE
`Atoms`. Requires the optional dependencies installed via
`pip install moyopy[interface]`. The magnetic counterpart
[`moyopy.interface.MoyoNonCollinearMagneticAdapter`][moyopy.interface.MoyoNonCollinearMagneticAdapter]
lives on the [Magnetic Space Group](api_magnetic_space_group.md) page.

::: moyopy.interface.MoyoAdapter
