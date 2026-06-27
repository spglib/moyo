# Examples

## Basic Usage

This example demonstrates the basic usage of the `moyopy` package.
First, we create a [`moyopy.Cell`][moyopy.Cell] representing a crystal structure, and then create a [`moyopy.MoyoDataset`][moyopy.MoyoDataset].
The [`moyopy.MoyoDataset`][moyopy.MoyoDataset] contains symmetry information of the input crystal structure: for example, the space group number, symmetry operations, and standardized cell.
When we need a secondary symmetry information such as Hermann-Mauguin symbol for the space group type, we can use the [`moyopy.HallSymbolEntry`][moyopy.HallSymbolEntry] to access the symmetry database.

```python
--8<-- "examples/basic.py"
```

## Accessing space-group type information

You can access the space-group classification information using [`moyopy.SpaceGroupType`][moyopy.SpaceGroupType].

```python
--8<-- "examples/space_group_type.py"
```

You can reverse-lookup an ITA space group number from a Hermann-Mauguin symbol by looping over [`moyopy.SpaceGroupType`][moyopy.SpaceGroupType].

```python
--8<-- "examples/hm_to_number.py"
```
