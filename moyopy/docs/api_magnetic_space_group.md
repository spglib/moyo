# Magnetic Space Group API

Magnetic crystal structures, symmetry analysis, classification tables, and
group identification for the 1,651 magnetic space groups. Shared types such
as [`moyopy.Cell`][moyopy.Cell], [`moyopy.Operations`][moyopy.Operations],
[`moyopy.UnimodularTransformation`][moyopy.UnimodularTransformation], and
[`moyopy.PointGroup`][moyopy.PointGroup] are documented on the
[API Reference](api.md) hub.

## Magnetic core types

Lattice + sites with magnetic moments, plus the magnetic operation
container.

::: moyopy.CollinearMagneticCell

::: moyopy.NonCollinearMagneticCell

::: moyopy.MagneticOperations

## Symmetry datasets

Run a symmetry analysis on a magnetic cell and inspect the result.

::: moyopy.MoyoCollinearMagneticDataset

::: moyopy.MoyoNonCollinearMagneticDataset

## Crystallographic data

Classification tables and helpers to fetch operations by UNI number.

::: moyopy.MagneticSpaceGroupType

::: moyopy.magnetic_operations_from_uni_number

## Group identification

Identify magnetic space groups from a primitive list of magnetic symmetry
operations.

::: moyopy.MagneticSpaceGroup

## Adapters

Convert between [`moyopy.NonCollinearMagneticCell`][moyopy.NonCollinearMagneticCell]
and pymatgen `Structure` / ASE `Atoms`. Requires the optional dependencies
installed via `pip install moyopy[interface]`.

::: moyopy.interface.MoyoNonCollinearMagneticAdapter
