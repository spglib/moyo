# Layer Group API

```{eval-rst}
.. currentmodule:: moyopy
```

Symmetry analysis and classification tables for the 80 layer groups. Shared
types such as {py:class}`moyopy.Cell`, {py:class}`moyopy.Operations`,
{py:class}`moyopy.UnimodularTransformation`, and {py:class}`moyopy.PointGroup`
are documented on the [API Reference](api.md) hub.

## Symmetry datasets

Run a symmetry analysis on a {py:class}`moyopy.Cell` treated as a layer (slab)
and inspect the result.

```{eval-rst}
.. autoapisummary::

   moyopy.MoyoLayerDataset

.. autoapiclass:: moyopy.MoyoLayerDataset
   :members:
```

## Crystallographic data

Layer-group settings, classification tables, and helpers to fetch operations
by layer-group number.

```{eval-rst}
.. autoapisummary::

   moyopy.LayerSetting
   moyopy.LayerCentering
   moyopy.LayerHallSymbolEntry
   moyopy.LayerGroupType
   moyopy.LayerArithmeticCrystalClass
   moyopy.operations_from_layer_number

.. autoapiclass:: moyopy.LayerSetting
   :members:

.. autoapiclass:: moyopy.LayerCentering
   :members:

.. autoapiclass:: moyopy.LayerHallSymbolEntry
   :members:

.. autoapiclass:: moyopy.LayerGroupType
   :members:

.. autoapiclass:: moyopy.LayerArithmeticCrystalClass
   :members:

.. autoapifunction:: moyopy.operations_from_layer_number
```
