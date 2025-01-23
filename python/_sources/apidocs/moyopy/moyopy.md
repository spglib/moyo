# {py:mod}`moyopy`

```{py:module} moyopy
```

```{autodoc2-docstring} moyopy
:parser: rst
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`MoyoDataset <moyopy.MoyoDataset>`
  - ```{autodoc2-docstring} moyopy.MoyoDataset
    :parser: rst
    :summary:
    ```
````

### Data

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`__version__ <moyopy.__version__>`
  - ```{autodoc2-docstring} moyopy.__version__
    :parser: rst
    :summary:
    ```
* - {py:obj}`__all__ <moyopy.__all__>`
  - ```{autodoc2-docstring} moyopy.__all__
    :parser: rst
    :summary:
    ```
````

### API

````{py:data} __version__
:canonical: moyopy.__version__
:type: str
:value: >
   None

```{autodoc2-docstring} moyopy.__version__
:parser: rst
```

````

`````{py:class} MoyoDataset(cell: moyopy._base.Cell, symprec: float = 0.0001, angle_tolerance: float | None = None, setting: moyopy._data.Setting | None = None)
:canonical: moyopy.MoyoDataset

```{autodoc2-docstring} moyopy.MoyoDataset
:parser: rst
```

```{rubric} Initialization
```

```{autodoc2-docstring} moyopy.MoyoDataset.__init__
:parser: rst
```

````{py:property} number
:canonical: moyopy.MoyoDataset.number
:type: int

```{autodoc2-docstring} moyopy.MoyoDataset.number
:parser: rst
```

````

````{py:property} hall_number
:canonical: moyopy.MoyoDataset.hall_number
:type: int

```{autodoc2-docstring} moyopy.MoyoDataset.hall_number
:parser: rst
```

````

````{py:property} operations
:canonical: moyopy.MoyoDataset.operations
:type: moyopy._base.Operations

```{autodoc2-docstring} moyopy.MoyoDataset.operations
:parser: rst
```

````

````{py:property} orbits
:canonical: moyopy.MoyoDataset.orbits
:type: list[int]

```{autodoc2-docstring} moyopy.MoyoDataset.orbits
:parser: rst
```

````

````{py:property} wyckoffs
:canonical: moyopy.MoyoDataset.wyckoffs
:type: list[str]

```{autodoc2-docstring} moyopy.MoyoDataset.wyckoffs
:parser: rst
```

````

````{py:property} site_symmetry_symbols
:canonical: moyopy.MoyoDataset.site_symmetry_symbols
:type: list[str]

```{autodoc2-docstring} moyopy.MoyoDataset.site_symmetry_symbols
:parser: rst
```

````

````{py:property} std_cell
:canonical: moyopy.MoyoDataset.std_cell
:type: moyopy._base.Cell

```{autodoc2-docstring} moyopy.MoyoDataset.std_cell
:parser: rst
```

````

````{py:property} std_linear
:canonical: moyopy.MoyoDataset.std_linear
:type: list[list[float]]

```{autodoc2-docstring} moyopy.MoyoDataset.std_linear
:parser: rst
```

````

````{py:property} std_origin_shift
:canonical: moyopy.MoyoDataset.std_origin_shift
:type: list[float]

```{autodoc2-docstring} moyopy.MoyoDataset.std_origin_shift
:parser: rst
```

````

````{py:property} std_rotation_matrix
:canonical: moyopy.MoyoDataset.std_rotation_matrix
:type: list[list[float]]

```{autodoc2-docstring} moyopy.MoyoDataset.std_rotation_matrix
:parser: rst
```

````

````{py:property} prim_std_cell
:canonical: moyopy.MoyoDataset.prim_std_cell
:type: moyopy._base.Cell

```{autodoc2-docstring} moyopy.MoyoDataset.prim_std_cell
:parser: rst
```

````

````{py:property} prim_std_linear
:canonical: moyopy.MoyoDataset.prim_std_linear
:type: list[list[float]]

```{autodoc2-docstring} moyopy.MoyoDataset.prim_std_linear
:parser: rst
```

````

````{py:property} prim_std_origin_shift
:canonical: moyopy.MoyoDataset.prim_std_origin_shift
:type: list[float]

```{autodoc2-docstring} moyopy.MoyoDataset.prim_std_origin_shift
:parser: rst
```

````

````{py:property} mapping_std_prim
:canonical: moyopy.MoyoDataset.mapping_std_prim
:type: list[int]

```{autodoc2-docstring} moyopy.MoyoDataset.mapping_std_prim
:parser: rst
```

````

````{py:property} symprec
:canonical: moyopy.MoyoDataset.symprec
:type: float

```{autodoc2-docstring} moyopy.MoyoDataset.symprec
:parser: rst
```

````

````{py:property} angle_tolerance
:canonical: moyopy.MoyoDataset.angle_tolerance
:type: float | None

```{autodoc2-docstring} moyopy.MoyoDataset.angle_tolerance
:parser: rst
```

````

`````

````{py:data} __all__
:canonical: moyopy.__all__
:value: >
   ['Cell', 'Operations', 'Setting', 'Centering', 'HallSymbolEntry', 'operations_from_number', '__versi...

```{autodoc2-docstring} moyopy.__all__
:parser: rst
```

````
