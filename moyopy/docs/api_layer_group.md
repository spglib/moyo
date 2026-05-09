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

Layer-group settings.

```{eval-rst}
.. autoapisummary::

   moyopy.LayerSetting

.. autoapiclass:: moyopy.LayerSetting
   :members:
```
