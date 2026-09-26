# Migration Guide from spglib

Based on spglib v2.6.0

## Symmetry tolerance

### `symprec` and `mag_symprec`

Same as in spglib.

### `angle_tolerance`

`angle_tolerance` in moyopy is in **radian** unit, while in spglib it is in **degree** unit.
That being said, moyopy is expected to work without specifying `angle_tolerance` in most cases.

## Space group

### Cell representation

Spglib's `cell = (lattice, positions, numbers)` corresponds to `moyopy.Cell(basis=lattice, positions=positions, numbers=numbers)`.

### Atomic equivalence

`MoyoDataset.orbits` corresponds to spglib's
[`crystallographic_orbits`](https://spglib.readthedocs.io/en/v2.7.0/dataset.html#crystallographic-orbits).
Both group the input atoms using the symmetry of the primitive cell. Spglib's
`equivalent_atoms` uses the symmetry operations found for the input cell.

These partitions usually agree, but the input-cell symmetry can be lower for
some supercells. In that case, `equivalent_atoms` may split a crystallographic
orbit into smaller groups. `MoyoDataset.orbits` is therefore not a general
replacement for `equivalent_atoms`, which moyopy does not directly expose.

When operations have non-integer rotation matrices in the input-cell basis,
moyo omits them and warns about this possible difference.

### Space-group symmetry search

#### `spglib.get_symmetry()`

Replace with `MoyoDataset`.

Let `symmetry` be the returned dictionary of `spglib.get_symmetry()`, the fields correspond to `MoyoDataset` fields as follows:

- `symmetry['rotations']` -> `MoyoDataset.operations.rotations`
- `symmetry['translations']` -> `MoyoDataset.operations.translations`
- `symmetry['equivalent_atoms']` is not directly available in `MoyoDataset`; see
  [Atomic equivalence](#atomic-equivalence) for the distinction from `MoyoDataset.orbits`.

#### `spglib.get_symmetry_dataset()`

Replace with `MoyoDataset`.

If you need to specify `hall_number`, use

```text
MoyoDataset(..., setting=Setting.hall_number(hall_number))
```

`spglib.SpglibDataset` fields correspond to `MoyoDataset` fields as follows:

- Space-group type
  - `number` -> `MoyoDataset.number`
  - `hall_number` -> `MoyoDataset.hall_number`
- Symmetry operations in the input cell
  - `rotations` -> `MoyoDataset.operations.rotations`
  - `translations` -> `MoyoDataset.operations.translations`
- Site symmetry
  - `crystallographic_orbits` -> `MoyoDataset.orbits`
  - `wyckoffs` -> `MoyoDataset.wyckoffs`
  - `site_symmetry_symbols` -> `MoyoDataset.site_symmetry_symbols`
- Standardized cell
  - `std_lattice` -> `MoyoDataset.std_cell.basis`
  - `std_positions` -> `MoyoDataset.std_cell.positions`
  - `std_types` -> `MoyoDataset.std_cell.numbers`
  - `transformation_matrix` -> `MoyoDataset.std_linear`
  - `origin_shift` -> `MoyoDataset.std_origin_shift`
  - `std_rotation_matrix` -> `MoyoDataset.std_rotation_matrix`
- Primitive standardized cell
  - `mapping_to_primitive` -> `MoyoDataset.mapping_std_prim`
- Secondary information
  - `hall` -> Obtain `hall_number` from `MoyoDataset.hall_number` and Access `HallSymbolEntry(hall_number).hm_short`
  - `choice` -> Obtain `hall_number` from `MoyoDataset.hall_number` and Access `HallSymbolEntry(hall_number).centering`
  - `international` -> Obtain `number` from `MoyoDataset.number` and Access `SpaceGroupType(number).hm_short`
- Not supported in moyopy
  - `equivalent_atoms`: See [Atomic equivalence](#atomic-equivalence) for the distinction from `MoyoDataset.orbits`
  - `primitive_lattice`: Use `MoyoDataset.prim_std_cell` if you need a primitive standardized cell
  - `std_mapping_to_primitive`: Recreate `MoyoDataset` again with `MoyoDataset.prim_std_cell`
  - `pointgroup`

### Space-group type search

#### `spglib.get_spacegroup()`

Moyopy does not directly support this function.
[See `spglib.get_symmetry_dataset()`](#spglibget_symmetry_dataset).

### Space-group dataset access

#### `spglib.get_symmetry_from_database()`

Replace `spglib.get_symmetry_from_database(hall_number)` with `operations_from_number(number, setting)`.
An ITA number, `number`, can be obtained from `HallSymbolEntry(hall_number).number`.

#### `spglib.get_spacegroup_type()`

Replace `spglib.get_spacegroup_type(hall_number)` with `SpaceGroupType(number)`.
An ITA number, `number`, can be obtained from `HallSymbolEntry(hall_number).number`.

`spglib.SpaceGroupType` fields correspond to `SpaceGroupType` fields as follows:

- Space-group type
  - `number` -> `SpaceGroupType.number`
  - `international_short` -> `SpaceGroupType.hm_short`
  - `international_full` -> `SpaceGroupType.hm_full`
  - `international` -> `SpaceGroupType.hm_short`
- Arithmetic crystal class
  - `arithmetic_crystal_class_number` -> `SpaceGroupType.arithmetic_number`
  - `arithmetic_crystal_class_symbol` -> `SpaceGroupType.arithmetic_symbol`
- Other classifications
  - `pointgroup_international` -> `SpaceGroupType.geometric_crystal_class`
- Not supported in moyopy
  - Hall symbol information (`choice`, `hall_number`, `hall_symbol`) -> Use `HallSymbolEntry` instead
  - `pointgroup_schoenflies`
  - `schoenflies`

#### `spglib.get_spacegroup_type_from_symmetry()`

Replace `spglib.get_spacegroup_type_from_symmetry(rotations, translations, lattice)` it with `SpaceGroup(prim_rotations, prim_translations, basis=lattice)`.
Be careful that `SpaceGroup` works only with symmetry operations in primitive.

## Magnetic symmetry

### Magnetic cell representation

Spglib's magnetic cell `cell = (lattice, positions, numbers, magmoms)` with collinear (N, )-shape `magmoms` corresponds to `MoyoCollinearMagneticCell(basis=lattice, positions=positions, numbers=numbers, magnetic_moments=magmoms)`.
Spglib's magnetic cell `cell = (lattice, positions, numbers, magmoms)` with noncollinear (N, 3)-shape `magmoms` corresponds to `MoyoNonCollinearMagneticCell(basis=lattice, positions=positions, numbers=numbers, magnetic_moments=magmoms)`.

### `spglib.get_magnetic_symmetry()`

Replace with `MoyoCollinearMagneticDataset` for collinear and `MoyoNonCollinearMagneticDataset` for noncollinear.

Let `symmetry` be the returned dictionary of `spglib.get_magnetic_symmetry()`, the fields correspond to `MoyoNonCollinearMagneticDataset` as follows:

- `symmetry['rotations']` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.rotations`
- `symmetry['translations']` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.translations`
- `symmetry['time_reversals']` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.time_reversals`
- `symmetry['equivalent_atoms']` is not directly available in `MoyoNonCollinearMagneticDataset`. However, `MoyoNonCollinearMagneticDataset.orbits` gives Spglib's `crystallographic_orbits`.
- `symmetry['primitive_lattice']` -> `MoyoNonCollinearMagneticDataset.prim_std_mag_cell.basis`
  Reread `MoyoNonCollinearMagneticDataset` as `MoyoCollinearMagneticDataset` for collinear case.

### `spglib.get_magnetic_symmetry_dataset()`

Replace with `MoyoCollinearMagneticDataset` for collinear and `MoyoNonCollinearMagneticDataset` for noncollinear.

`spglib.SpglibMagneticDataset` fields correspond to `MoyoNonCollinearMagneticDataset` fields as follows:

- Magnetic space-group type
  - `uni_number` -> `MoyoNonCollinearMagneticDataset.uni_number`
- Magnetic symmetry operations in the input cell
  - `rotations` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.rotations`
  - `translations` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.translations`
  - `time_reversals` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.time_reversals`
  - `n_operations` -> `MoyoNonCollinearMagneticDataset.magnetic_operations.num_operations`
- Site symmetry
  - `equivalent_atoms` -> Use `MoyoNonCollinearMagneticDataset.orbits` to get crystallographic orbits instead
- Standardized magnetic cell
  - `std_lattice` -> `MoyoNonCollinearMagneticDataset.std_mag_cell.basis`
  - `std_positions` -> `MoyoNonCollinearMagneticDataset.std_mag_cell.positions`
  - `std_types` -> `MoyoNonCollinearMagneticDataset.std_mag_cell.numbers`
  - `std_tensors` -> `MoyoNonCollinearMagneticDataset.std_mag_cell.magnetic_moments`
  - `n_std_atoms` -> `MoyoNonCollinearMagneticDataset.std_mag_cell.num_atoms`
  - `transformation_matrix` -> `MoyoNonCollinearMagneticDataset.std_linear`
  - `origin_shift` -> `MoyoNonCollinearMagneticDataset.std_origin_shift`
  - `std_rotation_matrix` -> `MoyoNonCollinearMagneticDataset.std_rotation_matrix`
- Secondary information
  - `msg_type` -> Obtain `uni_number` from `MoyoNonCollinearMagneticDataset.uni_number` and Access `MagneticSpaceGroupType(uni_number).construct_type`
- Not supported in moyopy
  - `primitive_lattice` -> Use `MoyoNonCollinearMagneticDataset.prim_std_mag_cell` if you need a primitive standardized magnetic cell
  - `hall_number` -> Obtain `uni_number` from `MoyoNonCollinearMagneticDataset.uni_number` and Access `MagneticSpaceGroupType(uni_number)` instead
  - `tensor_rank`
  - `n_atoms`
    Reread `MoyoNonCollinearMagneticDataset` as `MoyoCollinearMagneticDataset` for collinear case.

### `spglib.get_magnetic_spacegroup_type()`

Replace with `MagneticSpaceGroupType(uni_number)`.

`spglib.MagneticSpaceGroupType` fields correspond to `MagneticSpaceGroupType` fields as follows:

- `uni_number` -> `MagneticSpaceGroupType.uni_number`
- `litvin_number` -> `MagneticSpaceGroupType.litvin_number`
- `bns_number` -> `MagneticSpaceGroupType.bns_number`
- `og_number` -> `MagneticSpaceGroupType.og_number`
- `number` -> `MagneticSpaceGroupType.number`
- `type` -> `MagneticSpaceGroupType.construct_type`

### `spglib.get_magnetic_spacegroup_type_from_symmetry()`

Replace `spglib.get_magnetic_spacegroup_type_from_symmetry(rotations, translations, time_reversals, lattice)` with `MagneticSpaceGroup(prim_rotations, prim_translations, prim_time_reversals, basis=lattice)`.
Be careful that `MagneticSpaceGroup` works only with symmetry operations in primitive.

### `spglib.get_magnetic_symmetry_from_database()`

Replace it with `magnetic_operations_from_uni_number(uni_number)`.

## Standardization and finding primitive cell

### `spglib.standardize_cell()`

Create a `MoyoDataset` and read the desired standardized cell:

| Desired cell | spglib | moyopy |
| --- | --- | --- |
| Primitive standardized cell | `to_primitive=True` | `dataset.prim_std_cell` |
| Conventional standardized cell | `to_primitive=False` | `dataset.std_cell` |

To avoid the additional Cartesian rotation of the standardized lattice basis,
as with spglib's `no_idealize=True`, pass `rotate_basis=False` to `MoyoDataset`.
Use `rotate_basis=True` to apply the rotation returned as `std_rotation_matrix`.

Both libraries still standardize the cell when rotation is disabled. Moyo
refines atomic positions to match the detected symmetry for either value of
`rotate_basis`. If you need the original atomic positions, use your original
input cell.

Spglib's `no_idealize=True` skips its
[idealization step](https://spglib.readthedocs.io/en/v2.7.0/api.html#spg-standardize-cell),
but primitive-cell reduction can still
[average translation-equivalent positions](https://github.com/spglib/spglib/blob/v2.7.0/src/cell.c#L399-L460).
This can occur even with conventional output (`to_primitive=False`), so the flag
does not guarantee unchanged atomic positions.

### `spglib.find_primitive()`

Moyopy does not directly support this function.
[See `spglib.standardize_cell()`](#spglibstandardize_cell).

### `spglib.refine_cell()`

Moyopy does not directly support this function.
[See `spglib.standardize_cell()`](#spglibstandardize_cell).

## Lattice reduction

### `spglib.niggli_reduce()`

Use [`moyopy.niggli_reduce`][moyopy.niggli_reduce]. Both input and output bases
contain lattice vectors as rows. Moyopy also returns the integer transformation
matrix:

```python
import numpy as np
from moyopy import niggli_reduce

basis = np.array([[1.0, 0.0, 0.0], [5.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
reduced_basis, transformation = niggli_reduce(basis.tolist())
assert np.allclose(reduced_basis, np.array(transformation).T @ basis)
```

Moyopy uses the Rust implementation's fixed absolute tolerance and has no
`eps` argument. Choose units with typical lattice-vector lengths near one;
reduction and its predicates are not invariant under arbitrary rescaling.
Transformation coefficients must fit signed 32-bit integers; more extreme
shears are unsupported. Invalid input or a reported reduction failure raises
`ValueError`.
The choice among equivalent reduced bases can differ from spglib, particularly
on reduction boundaries; exact agreement with spglib's chosen cell is not
guaranteed.

### `spglib.delaunay_reduce()`

Use [`moyopy.delaunay_reduce`][moyopy.delaunay_reduce] with the same input and
return conventions:

```python
from moyopy import delaunay_reduce

reduced_basis, transformation = delaunay_reduce(basis.tolist())
```

As with Niggli reduction, there is no `eps` argument and invalid input raises
`ValueError`. Moyopy also provides [`moyopy.minkowski_reduce`][moyopy.minkowski_reduce]
and the predicates [`moyopy.is_niggli_reduced`][moyopy.is_niggli_reduced] and
[`moyopy.is_minkowski_reduced`][moyopy.is_minkowski_reduced].

## Kpoints: `spglib.get_ir_reciprocal_mesh()`

moyopy does not plan to support it
