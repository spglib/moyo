# API Reference

```{eval-rst}
.. currentmodule:: moyopy
```

The public API mirrors the Rust crate layout: `base` (crystal structures and
operation containers), `dataset` (symmetry analysis results), `data`
(crystallographic classification tables), and `identify` (group identification
primitives). The :mod:`moyopy.interface` adapters provide conversion helpers
to/from pymatgen and ASE.

## Crystal structures

Lattice + sites with optional magnetic moments, plus symmetry-operation
containers.

```{eval-rst}
.. autoapisummary::

   moyopy.Cell
   moyopy.CollinearMagneticCell
   moyopy.NonCollinearMagneticCell
   moyopy.Operations
   moyopy.MagneticOperations
   moyopy.UnimodularTransformation

.. autoapiclass:: moyopy.Cell
   :members:

.. autoapiclass:: moyopy.CollinearMagneticCell
   :members:

.. autoapiclass:: moyopy.NonCollinearMagneticCell
   :members:

.. autoapiclass:: moyopy.Operations
   :members:

.. autoapiclass:: moyopy.MagneticOperations
   :members:

.. autoapiclass:: moyopy.UnimodularTransformation
   :members:
```

## Symmetry datasets

Run a symmetry analysis on a :class:`Cell` or magnetic cell and inspect the
result.

```{eval-rst}
.. autoapisummary::

   moyopy.MoyoDataset
   moyopy.MoyoLayerDataset
   moyopy.MoyoCollinearMagneticDataset
   moyopy.MoyoNonCollinearMagneticDataset

.. autoapiclass:: moyopy.MoyoDataset
   :members:

.. autoapiclass:: moyopy.MoyoLayerDataset
   :members:

.. autoapiclass:: moyopy.MoyoCollinearMagneticDataset
   :members:

.. autoapiclass:: moyopy.MoyoNonCollinearMagneticDataset
   :members:
```

## Crystallographic data

Hall symbols, settings, centering, classification tables, and helpers to fetch
operations by ITA / UNI number.

```{eval-rst}
.. autoapisummary::

   moyopy.Setting
   moyopy.LayerSetting
   moyopy.Centering
   moyopy.HallSymbolEntry
   moyopy.SpaceGroupType
   moyopy.MagneticSpaceGroupType
   moyopy.ArithmeticCrystalClass
   moyopy.operations_from_number
   moyopy.magnetic_operations_from_uni_number

.. autoapiclass:: moyopy.Setting
   :members:

.. autoapiclass:: moyopy.LayerSetting
   :members:

.. autoapiclass:: moyopy.Centering
   :members:

.. autoapiclass:: moyopy.HallSymbolEntry
   :members:

.. autoapiclass:: moyopy.SpaceGroupType
   :members:

.. autoapiclass:: moyopy.MagneticSpaceGroupType
   :members:

.. autoapiclass:: moyopy.ArithmeticCrystalClass
   :members:

.. autoapifunction:: moyopy.operations_from_number

.. autoapifunction:: moyopy.magnetic_operations_from_uni_number
```

## Group identification

Identify point groups, space groups, and magnetic space groups from a primitive
list of symmetry operations.

```{eval-rst}
.. autoapisummary::

   moyopy.PointGroup
   moyopy.SpaceGroup
   moyopy.MagneticSpaceGroup
   moyopy.integral_normalizer

.. autoapiclass:: moyopy.PointGroup
   :members:

.. autoapiclass:: moyopy.SpaceGroup
   :members:

.. autoapiclass:: moyopy.MagneticSpaceGroup
   :members:

.. autoapifunction:: moyopy.integral_normalizer
```

## Adapters

Convert between :class:`moyopy.Cell` / :class:`moyopy.NonCollinearMagneticCell`
and pymatgen `Structure` / ASE `Atoms`. Requires the optional dependencies
installed via `pip install moyopy[interface]`.

```{eval-rst}
.. autoapisummary::

   moyopy.interface.MoyoAdapter
   moyopy.interface.MoyoNonCollinearMagneticAdapter

.. autoapiclass:: moyopy.interface.MoyoAdapter
   :members:

.. autoapiclass:: moyopy.interface.MoyoNonCollinearMagneticAdapter
   :members:
```
