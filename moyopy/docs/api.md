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

## Lattice reduction

Reduce a lattice given by three row-wise basis vectors. Each reduction returns
`(reduced_basis, transformation)` as nested lists, with an integer transformation
matrix satisfying `reduced_basis = transformation.T @ basis` in NumPy notation.

These functions use the Rust implementations' fixed absolute tolerances. Choose
units with typical lattice-vector lengths near one; changing the length scale
can change whether a basis is considered reduced. Transformation coefficients
must fit signed 32-bit integers. Extremely skewed bases requiring larger
coefficients are unsupported.

::: moyopy.niggli_reduce

::: moyopy.delaunay_reduce

::: moyopy.minkowski_reduce

::: moyopy.is_niggli_reduced

::: moyopy.is_minkowski_reduced

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

## Subgroup enumeration

Enumerate subgroups of a space group while preserving their embedding in the
parent group.

::: moyopy.enumerate_translationengleiche_subgroups

::: moyopy.TranslationengleicheSubgroup

::: moyopy.TranslationengleicheSubgroupConjugate

::: moyopy.TranslationengleicheSubgroupConjugacyClass

::: moyopy.enumerate_klassengleiche_subgroups

::: moyopy.enumerate_klassengleiche_subgroups_by_index

::: moyopy.KlassengleicheSubgroup

::: moyopy.KlassengleicheSubgroupConjugate

::: moyopy.KlassengleicheSubgroupConjugacyClass

## Adapters

Convert between [`moyopy.Cell`][moyopy.Cell] and pymatgen `Structure` / ASE
`Atoms`. Requires the optional dependencies installed via
`pip install moyopy[interface]`. The magnetic counterpart
[`moyopy.interface.MoyoNonCollinearMagneticAdapter`][moyopy.interface.MoyoNonCollinearMagneticAdapter]
lives on the [Magnetic Space Group](api_magnetic_space_group.md) page.

::: moyopy.interface.MoyoAdapter
