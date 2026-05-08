# Magnetic Space Group API

```{eval-rst}
.. currentmodule:: moyopy
```

Magnetic crystal structures, symmetry analysis, classification tables, and
group identification for the 1,651 magnetic space groups. Shared types such
as :class:`moyopy.Cell`, :class:`moyopy.Operations`,
:class:`moyopy.UnimodularTransformation`, and :class:`moyopy.PointGroup` are
documented on the [API Reference](api.md) hub.

## Magnetic core types

Lattice + sites with magnetic moments, plus the magnetic operation
container.

```{eval-rst}
.. autoapisummary::

   moyopy.CollinearMagneticCell
   moyopy.NonCollinearMagneticCell
   moyopy.MagneticOperations

.. autoapiclass:: moyopy.CollinearMagneticCell
   :members:

.. autoapiclass:: moyopy.NonCollinearMagneticCell
   :members:

.. autoapiclass:: moyopy.MagneticOperations
   :members:
```

## Symmetry datasets

Run a symmetry analysis on a magnetic cell and inspect the result.

```{eval-rst}
.. autoapisummary::

   moyopy.MoyoCollinearMagneticDataset
   moyopy.MoyoNonCollinearMagneticDataset

.. autoapiclass:: moyopy.MoyoCollinearMagneticDataset
   :members:

.. autoapiclass:: moyopy.MoyoNonCollinearMagneticDataset
   :members:
```

## Crystallographic data

Classification tables and helpers to fetch operations by UNI number.

```{eval-rst}
.. autoapisummary::

   moyopy.MagneticSpaceGroupType
   moyopy.magnetic_operations_from_uni_number

.. autoapiclass:: moyopy.MagneticSpaceGroupType
   :members:

.. autoapifunction:: moyopy.magnetic_operations_from_uni_number
```

## Group identification

Identify magnetic space groups from a primitive list of magnetic symmetry
operations.

```{eval-rst}
.. autoapisummary::

   moyopy.MagneticSpaceGroup

.. autoapiclass:: moyopy.MagneticSpaceGroup
   :members:
```

## Adapters

Convert between :class:`moyopy.NonCollinearMagneticCell` and pymatgen
`Structure` / ASE `Atoms`. Requires the optional dependencies installed via
`pip install moyopy[interface]`.

```{eval-rst}
.. autoapisummary::

   moyopy.interface.MoyoNonCollinearMagneticAdapter

.. autoapiclass:: moyopy.interface.MoyoNonCollinearMagneticAdapter
   :members:
```
