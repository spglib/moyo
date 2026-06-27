# Migration Guide from spglib

Based on spglib's C API v2.6.0.

This guide maps spglib's C functions (`spg_*`/`spgat_*`) and structs onto the moyoc
C API. For the moyopy (Python) equivalent, see the
[moyopy migration guide](https://spglib.github.io/moyo/python/migration/).

## General differences

### Basis vector layout

Spglib stores the lattice **column-wise**: `lattice[i][j]` is the `i`th Cartesian
component of the `j`th basis vector, so the `a` vector is
`(lattice[0][0], lattice[1][0], lattice[2][0])`.

moyoc stores the basis **row-wise**: `basis[i]` is the `i`th basis vector, so the
`a` vector is `(basis[0][0], basis[0][1], basis[0][2])`.

A moyoc basis is therefore the **transpose** of a spglib lattice (and vice versa
for every `*_lattice` output field):

```c
double basis[3][3];
for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        basis[i][j] = lattice[j][i]; // transpose
```

### Atomic numbers

Spglib's `types` (`int *`) corresponds to moyoc's `numbers` (`int32_t *`).

### Memory management

Spglib datasets are freed with `spg_free_dataset()`. Each moyoc object has its own
constructor/destructor pair (e.g. `moyo_dataset_new` / `moyo_dataset_free`). All
pointer fields inside an object are owned by that object and are released when you
free it; do not free them individually. Passing `NULL` to a `*_free` function is a
no-op, and every constructor returns `NULL` on failure.

## Symmetry tolerance

### `symprec` and `mag_symprec`

Same as in spglib.

### `angle_tolerance`

`angle_tolerance` in moyoc is in **radian** unit, while spglib's `spgat_*`
functions take it in **degree** unit. Pass a negative value to use the default
angle tolerance (moyoc is expected to work without specifying it in most cases).

### `hall_number`

Spglib selects a Hall setting through a `hall_number` argument. In moyoc, pass
`MOYO_SETTING_HALL_NUMBER` together with the desired `hall_number`:

```c
MoyoDataset *dataset = moyo_dataset_new(
    basis, positions, numbers, num_atoms,
    symprec, angle_tolerance, MOYO_SETTING_HALL_NUMBER, hall_number, rotate_basis
);
```

Otherwise use `MOYO_SETTING_SPGLIB` (smallest Hall number; spglib's default) or
`MOYO_SETTING_STANDARD`, and the `hall_number` argument is ignored.

## Space group

### `spg_get_symmetry()`

Replace with `moyo_dataset_new`, then read the symmetry operations from
`MoyoDataset.operations`:

- `rotation` -> `dataset->operations.rotations`
- `translation` -> `dataset->operations.translations`
- (number of operations) -> `dataset->operations.num_operations`
- `equivalent_atoms` is not directly available. Use `dataset->orbits`, which gives
  spglib's `crystallographic_orbits`.

### `spg_get_dataset()` and `spgat_get_dataset()`

Replace with `moyo_dataset_new` (free the result with `moyo_dataset_free`).

`SpglibDataset` fields correspond to `MoyoDataset` fields as follows:

- Space-group type
  - `spacegroup_number` -> `number`
  - `hall_number` -> `hall_number`
  - `international_symbol` -> `hm_symbol`
  - `pointgroup_symbol` -> Access `moyo_space_group_type_new(number)->geometric_crystal_class`
  - `hall_symbol` -> Access `moyo_hall_symbol_entry_new(hall_number)->hall_symbol`
  - `choice` -> Access `moyo_hall_symbol_entry_new(hall_number)->setting`
- Symmetry operations in the input cell
  - `rotations` -> `operations.rotations`
  - `translations` -> `operations.translations`
  - `n_operations` -> `operations.num_operations`
- Site symmetry
  - `crystallographic_orbits` -> `orbits`
  - `wyckoffs` -> `wyckoffs` (spglib returns integer indices; moyoc returns a
    NUL-terminated string of Wyckoff letters)
  - `site_symmetry_symbols` -> `site_symmetry_symbols`
  - `n_atoms` -> `num_atoms`
- Standardized cell
  - `std_lattice` -> `std_cell.basis` (transposed; see
    [Basis vector layout](#basis-vector-layout))
  - `std_positions` -> `std_cell.positions`
  - `std_types` -> `std_cell.numbers`
  - `n_std_atoms` -> `std_cell.num_atoms`
  - `transformation_matrix` -> `std_linear`
  - `origin_shift` -> `std_origin_shift`
  - `std_rotation_matrix` -> `std_rotation_matrix`
- Primitive standardized cell
  - `mapping_to_primitive` -> `mapping_std_prim`
- Not supported in moyoc
  - `equivalent_atoms` -> Use `orbits` to get crystallographic orbits instead
  - `primitive_lattice` -> Use `prim_std_cell` if you need a primitive
    standardized cell
  - `std_mapping_to_primitive` -> Recreate a `MoyoDataset` from `prim_std_cell`

### `spg_get_spacegroup()`

moyoc does not provide a single-call equivalent.
[See `spg_get_dataset()`](#spg_get_dataset-and-spgat_get_dataset) and read
`dataset->number` and `dataset->hm_symbol`.

### `spg_get_symmetry_from_database()`

Replace `spg_get_symmetry_from_database(rotations, translations, hall_number)` with
`moyo_operations_from_number(number, MOYO_SETTING_HALL_NUMBER, hall_number, primitive)`,
or pass an ITA `number` with another setting. An ITA number can be obtained from
`moyo_hall_symbol_entry_new(hall_number)->number`. Free the result with
`moyo_operations_free`.

### `spg_get_spacegroup_type()`

Replace `spg_get_spacegroup_type(hall_number)` with `moyo_space_group_type_new(number)`.
An ITA `number` can be obtained from `moyo_hall_symbol_entry_new(hall_number)->number`.

`SpglibSpacegroupType` fields correspond to `MoyoSpaceGroupType` fields as follows:

- Space-group type
  - `number` -> `number`
  - `international_short` -> `hm_short`
  - `international_full` -> `hm_full`
  - `international` -> `hm_short`
- Arithmetic crystal class
  - `arithmetic_crystal_class_number` -> `arithmetic_number`
  - `arithmetic_crystal_class_symbol` -> `arithmetic_symbol`
- Other classifications
  - `pointgroup_international` -> `geometric_crystal_class`
- Not supported in moyoc
  - Hall symbol information (`hall_number`, `hall_symbol`, `choice`) -> Use
    `moyo_hall_symbol_entry_new` instead
  - `schoenflies`
  - `pointgroup_schoenflies`

### `spg_get_spacegroup_type_from_symmetry()`

Replace
`spg_get_spacegroup_type_from_symmetry(rotations, translations, num_operations, lattice, symprec)`
with `moyo_space_group_new(prim_rotations, prim_translations, num_operations, basis, setting, hall_number, epsilon)`.
Be careful that `moyo_space_group_new` works only with symmetry operations in the
primitive cell. Free the result with `moyo_space_group_free`.

### `spg_get_pointgroup()`

Replace
`spg_get_pointgroup(symbol, trans_mat, rotations, num_rotations)` with
`moyo_point_group_new(prim_rotations, num_operations, basis)`, which works only with
rotations in the primitive cell.

- `trans_mat` -> `prim_trans_mat`
- The arithmetic crystal class number is `arithmetic_number`. Obtain the point-group
  symbol via `moyo_space_group_type_new(number)->geometric_crystal_class`.

Free the result with `moyo_point_group_free`.

## Magnetic symmetry

### Magnetic cell representation

Spglib passes magnetic moments through a separate `tensors` array. In moyoc:

- Collinear `(N,)` moments -> pass `magnetic_moments` (length `num_atoms`) to
  `moyo_collinear_magnetic_dataset_new`.
- Non-collinear `(N, 3)` moments -> pass `magnetic_moments` (Cartesian, length
  `3 * num_atoms`) to `moyo_noncollinear_magnetic_dataset_new`.

Spglib's `is_axial` corresponds to moyoc's `is_axial` argument.

### `spg_get_magnetic_symmetry()`

Replace with `moyo_collinear_magnetic_dataset_new` (collinear) or
`moyo_noncollinear_magnetic_dataset_new` (non-collinear), then read
`magnetic_operations`:

- `rotations` -> `magnetic_operations.rotations`
- `translations` -> `magnetic_operations.translations`
- `time_reversals` -> `magnetic_operations.time_reversals`
- (number of operations) -> `magnetic_operations.num_operations`
- `equivalent_atoms` is not directly available. Use `orbits`, which gives spglib's
  `crystallographic_orbits`.
- `primitive_lattice` -> `prim_std_mag_cell.basis` (transposed)

### `spg_get_magnetic_dataset()`

Replace with `moyo_collinear_magnetic_dataset_new` (collinear) or
`moyo_noncollinear_magnetic_dataset_new` (non-collinear). Free with the matching
`*_free` function.

`SpglibMagneticDataset` fields correspond to `MoyoNonCollinearMagneticDataset`
fields as follows (reread as `MoyoCollinearMagneticDataset` for the collinear case):

- Magnetic space-group type
  - `uni_number` -> `uni_number`
- Magnetic symmetry operations in the input cell
  - `rotations` -> `magnetic_operations.rotations`
  - `translations` -> `magnetic_operations.translations`
  - `time_reversals` -> `magnetic_operations.time_reversals`
  - `n_operations` -> `magnetic_operations.num_operations`
- Site symmetry
  - `equivalent_atoms` -> Use `orbits` to get crystallographic orbits instead
- Standardized magnetic cell
  - `std_lattice` -> `std_mag_cell.basis` (transposed)
  - `std_positions` -> `std_mag_cell.positions`
  - `std_types` -> `std_mag_cell.numbers`
  - `std_tensors` -> `std_mag_cell.magnetic_moments`
  - `n_std_atoms` -> `std_mag_cell.num_atoms`
  - `transformation_matrix` -> `std_linear`
  - `origin_shift` -> `std_origin_shift`
  - `std_rotation_matrix` -> `std_rotation_matrix`
- Secondary information
  - `msg_type` -> Access `moyo_magnetic_space_group_type_new(uni_number)->construct_type`
- Not supported in moyoc
  - `primitive_lattice` -> Use `prim_std_mag_cell` if you need a primitive
    standardized magnetic cell
  - `hall_number` -> Access `moyo_magnetic_space_group_type_new(uni_number)` instead
  - `tensor_rank`
  - `n_atoms`
  - `std_mapping_to_primitive`

### `spg_get_magnetic_spacegroup_type()`

Replace with `moyo_magnetic_space_group_type_new(uni_number)`.

`SpglibMagneticSpacegroupType` fields correspond to `MoyoMagneticSpaceGroupType`
fields as follows:

- `uni_number` -> `uni_number`
- `litvin_number` -> `litvin_number`
- `bns_number` -> `bns_number`
- `og_number` -> `og_number`
- `number` -> `number`
- `type` -> `construct_type`

### `spg_get_magnetic_spacegroup_type_from_symmetry()`

Replace
`spg_get_magnetic_spacegroup_type_from_symmetry(rotations, translations, time_reversals, num_operations, lattice, symprec)`
with
`moyo_magnetic_space_group_new(prim_rotations, prim_translations, prim_time_reversals, num_operations, basis, epsilon)`.
Be careful that `moyo_magnetic_space_group_new` works only with symmetry operations
in the primitive cell. Free with `moyo_magnetic_space_group_free`.

### `spg_get_magnetic_symmetry_from_database()`

Replace with `moyo_magnetic_operations_from_uni_number(uni_number, primitive)`. Free
the result with `moyo_magnetic_operations_free`.

## Standardization and finding primitive cell

### `spg_standardize_cell()`

`spg_standardize_cell(..., to_primitive, no_idealize, symprec)` is replaced with
creating a `MoyoDataset` and reading the appropriate cell. The `no_idealize` flag
is the **inverse** of moyoc's `rotate_basis` (idealizing the orientation
corresponds to `rotate_basis = true`):

- `to_primitive = 1` -> `moyo_dataset_new(..., !no_idealize)->prim_std_cell`
- `to_primitive = 0` -> `moyo_dataset_new(..., !no_idealize)->std_cell`

Remember to transpose `std_cell.basis` / `prim_std_cell.basis` back to spglib's
column-wise layout if you feed the result into spglib-style code.

### `spg_find_primitive()`

moyoc does not provide a single-call equivalent.
[See `spg_standardize_cell()`](#spg_standardize_cell).

### `spg_refine_cell()`

moyoc does not provide a single-call equivalent.
[See `spg_standardize_cell()`](#spg_standardize_cell).

## Lattice reduction

### `spg_niggli_reduce()`

Not supported yet.

### `spg_delaunay_reduce()`

Not supported yet.

## K-points: `spg_get_ir_reciprocal_mesh()`

moyoc does not plan to support it.
