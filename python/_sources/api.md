# API Reference

```{eval-rst}
.. currentmodule:: moyopy
```

The public API mirrors the Rust crate layout: `base` (crystal structures and
operation containers), `dataset` (symmetry analysis results), `data`
(crystallographic classification tables), and `identify` (group identification
primitives). The :mod:`moyopy.interface` adapters provide conversion helpers
to/from pymatgen and ASE.

This page documents the core types and the space-group API. Layer-group and
magnetic-space-group APIs live on additional pages:

```{toctree}
---
maxdepth: 1
---
Layer Group <api_layer_group>
Magnetic Space Group <api_magnetic_space_group>
```

## Core types

Lattice + sites, the non-magnetic operation container, and the unimodular
transformation type used across all group families.

```{eval-rst}
.. autoapisummary::

   moyopy.Cell
   moyopy.Operations
   moyopy.UnimodularTransformation

.. autoapiclass:: moyopy.Cell
   :members:

.. autoapiclass:: moyopy.Operations
   :members:

.. autoapiclass:: moyopy.UnimodularTransformation
   :members:
```

## Symmetry datasets

Run a symmetry analysis on a :class:`Cell` and inspect the result.

```{eval-rst}
.. autoapisummary::

   moyopy.MoyoDataset

.. autoapiclass:: moyopy.MoyoDataset
   :members:
```

## Crystallographic data

Hall symbols, settings, centering, classification tables, and helpers to fetch
operations by ITA number.

```{eval-rst}
.. autoapisummary::

   moyopy.Setting
   moyopy.Centering
   moyopy.HallSymbolEntry
   moyopy.SpaceGroupType
   moyopy.ArithmeticCrystalClass
   moyopy.operations_from_number

.. autoapiclass:: moyopy.Setting
   :members:

.. autoapiclass:: moyopy.Centering
   :members:

.. autoapiclass:: moyopy.HallSymbolEntry
   :members:

.. autoapiclass:: moyopy.SpaceGroupType
   :members:

.. autoapiclass:: moyopy.ArithmeticCrystalClass
   :members:

.. autoapifunction:: moyopy.operations_from_number
```

## Group identification

Identify point groups and space groups from a primitive list of symmetry
operations.

```{eval-rst}
.. autoapisummary::

   moyopy.PointGroup
   moyopy.SpaceGroup
   moyopy.integral_normalizer

.. autoapiclass:: moyopy.PointGroup
   :members:

.. autoapiclass:: moyopy.SpaceGroup
   :members:

.. autoapifunction:: moyopy.integral_normalizer
```

## Adapters

Convert between :class:`moyopy.Cell` and pymatgen `Structure` / ASE `Atoms`.
Requires the optional dependencies installed via
`pip install moyopy[interface]`. The magnetic counterpart
:class:`moyopy.interface.MoyoNonCollinearMagneticAdapter` lives on the
[Magnetic Space Group](api_magnetic_space_group.md) page.

```{eval-rst}
.. autoapisummary::

   moyopy.interface.MoyoAdapter

.. autoapiclass:: moyopy.interface.MoyoAdapter
   :members:
```
