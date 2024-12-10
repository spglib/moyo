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

* - {py:obj}`Cell <moyopy.Cell>`
  - ```{autodoc2-docstring} moyopy.Cell
    :parser: rst
    :summary:
    ```
* - {py:obj}`Operations <moyopy.Operations>`
  - ```{autodoc2-docstring} moyopy.Operations
    :parser: rst
    :summary:
    ```
* - {py:obj}`Setting <moyopy.Setting>`
  - ```{autodoc2-docstring} moyopy.Setting
    :parser: rst
    :summary:
    ```
* - {py:obj}`Centering <moyopy.Centering>`
  - ```{autodoc2-docstring} moyopy.Centering
    :parser: rst
    :summary:
    ```
* - {py:obj}`HallSymbolEntry <moyopy.HallSymbolEntry>`
  - ```{autodoc2-docstring} moyopy.HallSymbolEntry
    :parser: rst
    :summary:
    ```
* - {py:obj}`MoyoDataset <moyopy.MoyoDataset>`
  - ```{autodoc2-docstring} moyopy.MoyoDataset
    :parser: rst
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`operations_from_number <moyopy.operations_from_number>`
  - ```{autodoc2-docstring} moyopy.operations_from_number
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

`````{py:class} Cell(basis: list[list[float]], positions: list[list[float]], numbers: list[int])
:canonical: moyopy.Cell

```{autodoc2-docstring} moyopy.Cell
:parser: rst
```

```{rubric} Initialization
```

```{autodoc2-docstring} moyopy.Cell.__init__
:parser: rst
```

````{py:property} basis
:canonical: moyopy.Cell.basis
:type: list[list[float]]

```{autodoc2-docstring} moyopy.Cell.basis
:parser: rst
```

````

````{py:property} positions
:canonical: moyopy.Cell.positions
:type: list[list[float]]

```{autodoc2-docstring} moyopy.Cell.positions
:parser: rst
```

````

````{py:property} numbers
:canonical: moyopy.Cell.numbers
:type: list[int]

```{autodoc2-docstring} moyopy.Cell.numbers
:parser: rst
```

````

````{py:method} serialize_json() -> str
:canonical: moyopy.Cell.serialize_json

```{autodoc2-docstring} moyopy.Cell.serialize_json
:parser: rst
```

````

````{py:method} deserialize_json(json_str: str) -> moyopy.Cell
:canonical: moyopy.Cell.deserialize_json
:classmethod:

```{autodoc2-docstring} moyopy.Cell.deserialize_json
:parser: rst
```

````

`````

`````{py:class} Operations
:canonical: moyopy.Operations

```{autodoc2-docstring} moyopy.Operations
:parser: rst
```

````{py:property} rotations
:canonical: moyopy.Operations.rotations
:type: list[list[list[float]]]

```{autodoc2-docstring} moyopy.Operations.rotations
:parser: rst
```

````

````{py:property} translations
:canonical: moyopy.Operations.translations
:type: list[list[float]]

```{autodoc2-docstring} moyopy.Operations.translations
:parser: rst
```

````

````{py:property} num_operations
:canonical: moyopy.Operations.num_operations
:type: int

```{autodoc2-docstring} moyopy.Operations.num_operations
:parser: rst
```

````

````{py:method} __len__() -> int
:canonical: moyopy.Operations.__len__

```{autodoc2-docstring} moyopy.Operations.__len__
:parser: rst
```

````

`````

`````{py:class} Setting
:canonical: moyopy.Setting

```{autodoc2-docstring} moyopy.Setting
:parser: rst
```

````{py:method} spglib() -> moyopy.Setting
:canonical: moyopy.Setting.spglib
:classmethod:

```{autodoc2-docstring} moyopy.Setting.spglib
:parser: rst
```

````

````{py:method} standard() -> moyopy.Setting
:canonical: moyopy.Setting.standard
:classmethod:

```{autodoc2-docstring} moyopy.Setting.standard
:parser: rst
```

````

````{py:method} hall_number(hall_number: int) -> moyopy.Setting
:canonical: moyopy.Setting.hall_number
:classmethod:

```{autodoc2-docstring} moyopy.Setting.hall_number
:parser: rst
```

````

`````

````{py:function} operations_from_number(number: int, setting: moyopy.Setting) -> moyopy.Operations
:canonical: moyopy.operations_from_number

```{autodoc2-docstring} moyopy.operations_from_number
:parser: rst
```
````

````{py:class} Centering
:canonical: moyopy.Centering

```{autodoc2-docstring} moyopy.Centering
:parser: rst
```

````

`````{py:class} HallSymbolEntry(hall_number: int)
:canonical: moyopy.HallSymbolEntry

```{autodoc2-docstring} moyopy.HallSymbolEntry
:parser: rst
```

```{rubric} Initialization
```

```{autodoc2-docstring} moyopy.HallSymbolEntry.__init__
:parser: rst
```

````{py:property} hall_number
:canonical: moyopy.HallSymbolEntry.hall_number
:type: int

```{autodoc2-docstring} moyopy.HallSymbolEntry.hall_number
:parser: rst
```

````

````{py:property} number
:canonical: moyopy.HallSymbolEntry.number
:type: int

```{autodoc2-docstring} moyopy.HallSymbolEntry.number
:parser: rst
```

````

````{py:property} arithmetic_number
:canonical: moyopy.HallSymbolEntry.arithmetic_number
:type: int

```{autodoc2-docstring} moyopy.HallSymbolEntry.arithmetic_number
:parser: rst
```

````

````{py:property} setting
:canonical: moyopy.HallSymbolEntry.setting
:type: moyopy.Setting

```{autodoc2-docstring} moyopy.HallSymbolEntry.setting
:parser: rst
```

````

````{py:property} hall_symbol
:canonical: moyopy.HallSymbolEntry.hall_symbol
:type: str

```{autodoc2-docstring} moyopy.HallSymbolEntry.hall_symbol
:parser: rst
```

````

````{py:property} hm_short
:canonical: moyopy.HallSymbolEntry.hm_short
:type: str

```{autodoc2-docstring} moyopy.HallSymbolEntry.hm_short
:parser: rst
```

````

````{py:property} hm_full
:canonical: moyopy.HallSymbolEntry.hm_full
:type: str

```{autodoc2-docstring} moyopy.HallSymbolEntry.hm_full
:parser: rst
```

````

````{py:property} centering
:canonical: moyopy.HallSymbolEntry.centering
:type: moyopy.Centering

```{autodoc2-docstring} moyopy.HallSymbolEntry.centering
:parser: rst
```

````

`````

`````{py:class} MoyoDataset(cell: moyopy.Cell, symprec: float = 0.0001, angle_tolerance: float | None = None, setting: moyopy.Setting | None = None)
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
:type: moyopy.Operations

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
:type: moyopy.Cell

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
:type: moyopy.Cell

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
