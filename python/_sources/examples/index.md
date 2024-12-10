# Examples

## Basic Usage

This example demonstrates the basic usage of the `moyopy` package.
First, we create a {py:class}`moyopy.Cell` representing a crystal structure, and then create a {py:class}`moyopy.MoyoDataset`.
The {py:class}`moyopy.MoyoDataset` contains symmetry information of the input crystal structure: for example, the space group number, symmetry operations, and standardized cell.
When we need a secondary symmetry information such as Hermann-Mauguin symbol for the space group type, we can use the {py:class}`moyopy.HallSymbolEntry` to access the symmetry database.

```{literalinclude} ../../examples/basic.py
```
