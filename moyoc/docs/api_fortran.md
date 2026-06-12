# moyoc Fortran API reference

Fortran interface of [moyo](https://github.com/spglib/moyo) wrapping the
complete C API. See [the C API reference](./api.md) for the underlying
functions and types.

## Building

Configure moyoc with `-DMOYOC_BUILD_FORTRAN=ON` to build the Fortran interface
and expose the `moyo::moyof` CMake target. This requires a Fortran compiler
(for example gfortran) in addition to the usual moyoc prerequisites. Link
`moyo::moyof` into your target and `use moyo` from your Fortran sources; the
target carries the generated `moyo.mod` include directory.

```cmake
target_link_libraries(example PRIVATE moyo::moyof)
```

```fortran
use moyo
```

## Conventions

- Arrays cross the C boundary column-major. `basis(:, i)` is the `i`th basis
  vector and `positions(:, i)` are the fractional coordinates of the `i`th
  site, i.e. a transpose view of the row-major C arrays. The 3x3 matrices
  stored in datasets (`std_linear`, `rotations(:, :, k)`, and so on) are
  likewise transposed with respect to their C/Rust definition.
- Constructors return a `type(c_ptr)`. On failure they return NULL, which you
  test with `c_associated` (false means failure).
- Map a non-NULL handle to its matching `bind(c)` derived type with
  `c_f_pointer`, then read the fields of the resulting pointer.
- `c_ptr` components that point to arrays (for example `rotations`,
  `translations`, `orbits`, `mapping_std_prim`) are themselves mapped with
  `c_f_pointer` using the known length, for example
  `c_f_pointer(ops%rotations, rot, [3, 3, int(ops%num_operations)])`.
- Release each handle with the corresponding `*_free` routine. Do not free
  pointers embedded inside a dataset; free only the top-level handle that the
  constructor returned.
- Convert NUL-terminated C strings (the `c_ptr` string components) to Fortran
  `character` values with `moyo_to_string`.

## Setting constants

The module exposes the values of the C `MoyoSetting` and `MoyoLayerSetting`
enums as public integer parameters:

- `MOYO_SETTING_HALL_NUMBER = 0` -- Hall number specified by the `hall_number`
  argument.
- `MOYO_SETTING_SPGLIB = 1` -- the setting of the smallest Hall number.
- `MOYO_SETTING_STANDARD = 2` -- unique axis b, cell choice 1 for monoclinic,
  hexagonal axes for rhombohedral, and origin choice 2 for centrosymmetric
  space groups.
- `MOYO_LAYER_SETTING_HALL_NUMBER = 0` -- layer Hall number specified by the
  `hall_number` argument (1 - 116).
- `MOYO_LAYER_SETTING_SPGLIB = 1` -- the setting of the smallest layer Hall
  number for each layer group.
- `MOYO_LAYER_SETTING_STANDARD = 2` -- BCS standard choice (origin choice 2 for
  layer groups 52, 62, 64).

## Version

```fortran
function moyo_version() bind(c, name="moyo_version") result(version)
    type(c_ptr) :: version
end function moyo_version
```

Returns the moyoc version string as a `type(c_ptr)`; convert it with
`moyo_to_string`. The returned string is statically allocated and must not be
freed.

## Space-group dataset

```fortran
function moyo_dataset_new(basis, positions, numbers, num_atoms, symprec, &
                          angle_tolerance, setting, hall_number, rotate_basis) &
    result(dataset)
    real(c_double), intent(in) :: basis(3, 3)
    real(c_double), intent(in) :: positions(3, *)
    integer(c_int32_t), intent(in) :: numbers(*)
    integer(c_int32_t), value :: num_atoms
    real(c_double), value :: symprec
    real(c_double), value :: angle_tolerance
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    logical(c_bool), value :: rotate_basis
    type(c_ptr) :: dataset
end function moyo_dataset_new
```

Analyze the symmetry of the given cell. `basis(:, i)` is the `i`th basis
vector, `positions(:, i)` the fractional coordinates of site `i`, and
`numbers(i)` its atomic number, for `num_atoms` sites. `symprec` is the
distance tolerance; pass a negative `angle_tolerance` to use the default.
`setting` is one of the `MOYO_SETTING_*` constants and `hall_number` is only
used with `MOYO_SETTING_HALL_NUMBER`. Returns a pointer to a `MoyoDataset`, or
NULL on failure. Free the result with `moyo_dataset_free`.

```fortran
subroutine moyo_dataset_free(dataset)
    type(c_ptr), value :: dataset
end subroutine moyo_dataset_free
```

Free a dataset created by `moyo_dataset_new`.

### type `MoyoDataset`

- `number` (`integer(c_int32_t)`) -- ITA space-group number.
- `hall_number` (`integer(c_int32_t)`) -- Hall number.
- `hm_symbol` (`type(c_ptr)`) -- Hermann-Mauguin symbol C string; use
  `moyo_to_string`.
- `num_atoms` (`integer(c_int32_t)`) -- number of atoms in the input cell.
- `operations` (`type(MoyoOperations)`) -- symmetry operations in the input
  cell.
- `orbits` (`type(c_ptr)`) -- points to `num_atoms` int32 orbit indices.
- `wyckoffs` (`type(c_ptr)`) -- Wyckoff letters C string; use `moyo_to_string`.
- `site_symmetry_symbols` (`type(c_ptr)`) -- points to `num_atoms` C string
  pointers (`type(c_ptr)` array); convert each with `moyo_to_string`.
- `std_cell` (`type(MoyoCell)`) -- standardized conventional cell.
- `std_linear` (`real(c_double)(3, 3)`) -- transformation linear part
  (transpose view of the C matrix).
- `std_origin_shift` (`real(c_double)(3)`) -- origin shift of the
  standardization.
- `std_rotation_matrix` (`real(c_double)(3, 3)`) -- rigid rotation to the
  standardized cell (transpose view).
- `pearson_symbol` (`type(c_ptr)`) -- Pearson symbol C string; use
  `moyo_to_string`.
- `prim_std_cell` (`type(MoyoCell)`) -- standardized primitive cell.
- `prim_std_linear` (`real(c_double)(3, 3)`) -- primitive transformation linear
  part (transpose view).
- `prim_std_origin_shift` (`real(c_double)(3)`) -- primitive origin shift.
- `mapping_std_prim` (`type(c_ptr)`) -- points to `num_atoms` int32 indices
  mapping input sites into the standardized primitive cell.
- `symprec` (`real(c_double)`) -- distance tolerance used.
- `angle_tolerance` (`real(c_double)`) -- angle tolerance used.

## Layer-group dataset

```fortran
function moyo_layer_dataset_new(basis, positions, numbers, num_atoms, symprec, &
                                angle_tolerance, setting, hall_number, rotate_basis) &
    result(dataset)
    real(c_double), intent(in) :: basis(3, 3)
    real(c_double), intent(in) :: positions(3, *)
    integer(c_int32_t), intent(in) :: numbers(*)
    integer(c_int32_t), value :: num_atoms
    real(c_double), value :: symprec
    real(c_double), value :: angle_tolerance
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    logical(c_bool), value :: rotate_basis
    type(c_ptr) :: dataset
end function moyo_layer_dataset_new
```

Analyze the layer-group symmetry of the given cell. The third basis vector must
be the aperiodic stacking direction. `setting` is one of the
`MOYO_LAYER_SETTING_*` constants and `hall_number` is only used with
`MOYO_LAYER_SETTING_HALL_NUMBER`. Returns a pointer to a `MoyoLayerDataset`, or
NULL on failure. Free the result with `moyo_layer_dataset_free`.

```fortran
subroutine moyo_layer_dataset_free(dataset)
    type(c_ptr), value :: dataset
end subroutine moyo_layer_dataset_free
```

Free a dataset created by `moyo_layer_dataset_new`.

### type `MoyoLayerDataset`

The components mirror `MoyoDataset` exactly:

- `number` (`integer(c_int32_t)`) -- layer-group number.
- `hall_number` (`integer(c_int32_t)`) -- layer Hall number.
- `hm_symbol` (`type(c_ptr)`) -- Hermann-Mauguin symbol C string.
- `num_atoms` (`integer(c_int32_t)`) -- number of atoms in the input cell.
- `operations` (`type(MoyoOperations)`) -- symmetry operations in the input
  cell.
- `orbits` (`type(c_ptr)`) -- points to `num_atoms` int32 orbit indices.
- `wyckoffs` (`type(c_ptr)`) -- Wyckoff letters C string.
- `site_symmetry_symbols` (`type(c_ptr)`) -- points to `num_atoms` C string
  pointers.
- `std_cell` (`type(MoyoCell)`) -- standardized conventional cell.
- `std_linear` (`real(c_double)(3, 3)`) -- transformation linear part
  (transpose view).
- `std_origin_shift` (`real(c_double)(3)`) -- origin shift.
- `std_rotation_matrix` (`real(c_double)(3, 3)`) -- rigid rotation (transpose
  view).
- `pearson_symbol` (`type(c_ptr)`) -- Pearson symbol C string.
- `prim_std_cell` (`type(MoyoCell)`) -- standardized primitive cell.
- `prim_std_linear` (`real(c_double)(3, 3)`) -- primitive transformation linear
  part (transpose view).
- `prim_std_origin_shift` (`real(c_double)(3)`) -- primitive origin shift.
- `mapping_std_prim` (`type(c_ptr)`) -- points to `num_atoms` int32 mapping
  indices.
- `symprec` (`real(c_double)`) -- distance tolerance used.
- `angle_tolerance` (`real(c_double)`) -- angle tolerance used.

## Magnetic space-group datasets

```fortran
function moyo_collinear_magnetic_dataset_new(basis, positions, numbers, &
                                             magnetic_moments, num_atoms, symprec, &
                                             angle_tolerance, mag_symprec, is_axial, &
                                             rotate_basis) &
    result(dataset)
    real(c_double), intent(in) :: basis(3, 3)
    real(c_double), intent(in) :: positions(3, *)
    integer(c_int32_t), intent(in) :: numbers(*)
    real(c_double), intent(in) :: magnetic_moments(*)
    integer(c_int32_t), value :: num_atoms
    real(c_double), value :: symprec
    real(c_double), value :: angle_tolerance
    real(c_double), value :: mag_symprec
    logical(c_bool), value :: is_axial
    logical(c_bool), value :: rotate_basis
    type(c_ptr) :: dataset
end function moyo_collinear_magnetic_dataset_new
```

Analyze the magnetic symmetry of a collinear magnetic cell. `magnetic_moments`
holds one scalar moment per site (`num_atoms` values). Pass a negative
`mag_symprec` to reuse `symprec`. Returns a pointer to a
`MoyoCollinearMagneticDataset`, or NULL on failure. Free the result with
`moyo_collinear_magnetic_dataset_free`.

```fortran
subroutine moyo_collinear_magnetic_dataset_free(dataset)
    type(c_ptr), value :: dataset
end subroutine moyo_collinear_magnetic_dataset_free
```

Free a dataset created by `moyo_collinear_magnetic_dataset_new`.

```fortran
function moyo_noncollinear_magnetic_dataset_new(basis, positions, numbers, &
                                                magnetic_moments, num_atoms, symprec, &
                                                angle_tolerance, mag_symprec, is_axial, &
                                                rotate_basis) &
    result(dataset)
    real(c_double), intent(in) :: basis(3, 3)
    real(c_double), intent(in) :: positions(3, *)
    integer(c_int32_t), intent(in) :: numbers(*)
    real(c_double), intent(in) :: magnetic_moments(3, *)
    integer(c_int32_t), value :: num_atoms
    real(c_double), value :: symprec
    real(c_double), value :: angle_tolerance
    real(c_double), value :: mag_symprec
    logical(c_bool), value :: is_axial
    logical(c_bool), value :: rotate_basis
    type(c_ptr) :: dataset
end function moyo_noncollinear_magnetic_dataset_new
```

Analyze the magnetic symmetry of a non-collinear magnetic cell, where
`magnetic_moments(:, i)` is the Cartesian moment of site `i`. Pass a negative
`mag_symprec` to reuse `symprec`. Returns a pointer to a
`MoyoNonCollinearMagneticDataset`, or NULL on failure. Free the result with
`moyo_noncollinear_magnetic_dataset_free`.

```fortran
subroutine moyo_noncollinear_magnetic_dataset_free(dataset)
    type(c_ptr), value :: dataset
end subroutine moyo_noncollinear_magnetic_dataset_free
```

Free a dataset created by `moyo_noncollinear_magnetic_dataset_new`.

### type `MoyoCollinearMagneticDataset`

- `uni_number` (`integer(c_int32_t)`) -- UNI number of the magnetic space
  group.
- `num_atoms` (`integer(c_int32_t)`) -- number of atoms in the input cell.
- `magnetic_operations` (`type(MoyoMagneticOperations)`) -- magnetic symmetry
  operations.
- `orbits` (`type(c_ptr)`) -- points to `num_atoms` int32 orbit indices.
- `std_mag_cell` (`type(MoyoCollinearMagneticCell)`) -- standardized
  conventional magnetic cell.
- `std_linear` (`real(c_double)(3, 3)`) -- transformation linear part
  (transpose view).
- `std_origin_shift` (`real(c_double)(3)`) -- origin shift.
- `std_rotation_matrix` (`real(c_double)(3, 3)`) -- rigid rotation (transpose
  view).
- `prim_std_mag_cell` (`type(MoyoCollinearMagneticCell)`) -- standardized
  primitive magnetic cell.
- `prim_std_linear` (`real(c_double)(3, 3)`) -- primitive transformation linear
  part (transpose view).
- `prim_std_origin_shift` (`real(c_double)(3)`) -- primitive origin shift.
- `mapping_std_prim` (`type(c_ptr)`) -- points to `num_atoms` int32 mapping
  indices.
- `symprec` (`real(c_double)`) -- distance tolerance used.
- `angle_tolerance` (`real(c_double)`) -- angle tolerance used.
- `mag_symprec` (`real(c_double)`) -- magnetic moment tolerance used.

### type `MoyoNonCollinearMagneticDataset`

The components match `MoyoCollinearMagneticDataset` except that the cell fields
use the non-collinear cell type:

- `uni_number` (`integer(c_int32_t)`) -- UNI number.
- `num_atoms` (`integer(c_int32_t)`) -- number of atoms in the input cell.
- `magnetic_operations` (`type(MoyoMagneticOperations)`) -- magnetic symmetry
  operations.
- `orbits` (`type(c_ptr)`) -- points to `num_atoms` int32 orbit indices.
- `std_mag_cell` (`type(MoyoNonCollinearMagneticCell)`) -- standardized
  conventional magnetic cell.
- `std_linear` (`real(c_double)(3, 3)`) -- transformation linear part
  (transpose view).
- `std_origin_shift` (`real(c_double)(3)`) -- origin shift.
- `std_rotation_matrix` (`real(c_double)(3, 3)`) -- rigid rotation (transpose
  view).
- `prim_std_mag_cell` (`type(MoyoNonCollinearMagneticCell)`) -- standardized
  primitive magnetic cell.
- `prim_std_linear` (`real(c_double)(3, 3)`) -- primitive transformation linear
  part (transpose view).
- `prim_std_origin_shift` (`real(c_double)(3)`) -- primitive origin shift.
- `mapping_std_prim` (`type(c_ptr)`) -- points to `num_atoms` int32 mapping
  indices.
- `symprec` (`real(c_double)`) -- distance tolerance used.
- `angle_tolerance` (`real(c_double)`) -- angle tolerance used.
- `mag_symprec` (`real(c_double)`) -- magnetic moment tolerance used.

## Symmetry operations from the database

```fortran
function moyo_operations_from_number(number, setting, hall_number, primitive) &
    result(operations)
    integer(c_int32_t), value :: number
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    logical(c_bool), value :: primitive
    type(c_ptr) :: operations
end function moyo_operations_from_number
```

Symmetry operations for a space-group ITA number (1 - 230). `setting` is one of
the `MOYO_SETTING_*` constants. Returns a pointer to `MoyoOperations`, or NULL
on invalid input. Free the result with `moyo_operations_free`.

```fortran
function moyo_operations_from_layer_number(number, setting, hall_number, primitive) &
    result(operations)
    integer(c_int32_t), value :: number
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    logical(c_bool), value :: primitive
    type(c_ptr) :: operations
end function moyo_operations_from_layer_number
```

Symmetry operations for a layer-group number (1 - 80). `setting` is one of the
`MOYO_LAYER_SETTING_*` constants. Returns a pointer to `MoyoOperations`, or NULL
on invalid input. Free the result with `moyo_operations_free`.

```fortran
subroutine moyo_operations_free(operations)
    type(c_ptr), value :: operations
end subroutine moyo_operations_free
```

Free operations returned by `moyo_operations_from_number` or
`moyo_operations_from_layer_number`. Do not call this on operations embedded in
a dataset.

```fortran
function moyo_magnetic_operations_from_uni_number(uni_number, primitive) &
    result(magnetic_operations)
    integer(c_int32_t), value :: uni_number
    logical(c_bool), value :: primitive
    type(c_ptr) :: magnetic_operations
end function moyo_magnetic_operations_from_uni_number
```

Magnetic symmetry operations for a UNI number (1 - 1651). Returns a pointer to
`MoyoMagneticOperations`, or NULL on invalid input. Free the result with
`moyo_magnetic_operations_free`.

```fortran
subroutine moyo_magnetic_operations_free(magnetic_operations)
    type(c_ptr), value :: magnetic_operations
end subroutine moyo_magnetic_operations_free
```

Free magnetic operations returned by
`moyo_magnetic_operations_from_uni_number`. Do not call this on operations
embedded in a dataset.

### type `MoyoOperations`

- `rotations` (`type(c_ptr)`) -- points to `num_operations` int32 3x3 matrices;
  map with `c_f_pointer(ops%rotations, rot, [3, 3, int(ops%num_operations)])`,
  where each `rot(:, :, k)` is the transposed matrix.
- `translations` (`type(c_ptr)`) -- points to `num_operations` real translation
  vectors; map with shape `[3, int(ops%num_operations)]`.
- `num_operations` (`integer(c_int32_t)`) -- number of operations.

### type `MoyoMagneticOperations`

- `rotations` (`type(c_ptr)`) -- points to `num_operations` int32 3x3 matrices
  (as in `MoyoOperations`).
- `translations` (`type(c_ptr)`) -- points to `num_operations` real translation
  vectors.
- `time_reversals` (`type(c_ptr)`) -- points to `num_operations`
  `logical(c_bool)` time-reversal flags.
- `num_operations` (`integer(c_int32_t)`) -- number of operations.

## Group-type metadata

```fortran
function moyo_hall_symbol_entry_new(hall_number) result(entry)
    integer(c_int32_t), value :: hall_number
    type(c_ptr) :: entry
end function moyo_hall_symbol_entry_new
```

Hall symbol entry for `hall_number` (1 - 530). Returns a pointer to
`MoyoHallSymbolEntry`, or NULL if out of range. Free with
`moyo_hall_symbol_entry_free`.

```fortran
subroutine moyo_hall_symbol_entry_free(entry)
    type(c_ptr), value :: entry
end subroutine moyo_hall_symbol_entry_free
```

Free an entry created by `moyo_hall_symbol_entry_new`.

```fortran
function moyo_layer_hall_symbol_entry_new(hall_number) result(entry)
    integer(c_int32_t), value :: hall_number
    type(c_ptr) :: entry
end function moyo_layer_hall_symbol_entry_new
```

Layer Hall symbol entry for `hall_number` (1 - 116). Returns a pointer to
`MoyoLayerHallSymbolEntry`, or NULL if out of range. Free with
`moyo_layer_hall_symbol_entry_free`.

```fortran
subroutine moyo_layer_hall_symbol_entry_free(entry)
    type(c_ptr), value :: entry
end subroutine moyo_layer_hall_symbol_entry_free
```

Free an entry created by `moyo_layer_hall_symbol_entry_new`.

```fortran
function moyo_space_group_type_new(number) result(space_group_type)
    integer(c_int32_t), value :: number
    type(c_ptr) :: space_group_type
end function moyo_space_group_type_new
```

Space-group type information for ITA `number` (1 - 230). Returns a pointer to
`MoyoSpaceGroupType`, or NULL if out of range. Free with
`moyo_space_group_type_free`.

```fortran
subroutine moyo_space_group_type_free(space_group_type)
    type(c_ptr), value :: space_group_type
end subroutine moyo_space_group_type_free
```

Free a space-group type created by `moyo_space_group_type_new`.

```fortran
function moyo_layer_group_type_new(number) result(layer_group_type)
    integer(c_int32_t), value :: number
    type(c_ptr) :: layer_group_type
end function moyo_layer_group_type_new
```

Layer-group type information for layer-group `number` (1 - 80). Returns a
pointer to `MoyoLayerGroupType`, or NULL if out of range. Free with
`moyo_layer_group_type_free`.

```fortran
subroutine moyo_layer_group_type_free(layer_group_type)
    type(c_ptr), value :: layer_group_type
end subroutine moyo_layer_group_type_free
```

Free a layer-group type created by `moyo_layer_group_type_new`.

```fortran
function moyo_magnetic_space_group_type_new(uni_number) result(magnetic_space_group_type)
    integer(c_int32_t), value :: uni_number
    type(c_ptr) :: magnetic_space_group_type
end function moyo_magnetic_space_group_type_new
```

Magnetic space-group type information for `uni_number` (1 - 1651). Returns a
pointer to `MoyoMagneticSpaceGroupType`, or NULL if out of range. Free with
`moyo_magnetic_space_group_type_free`.

```fortran
subroutine moyo_magnetic_space_group_type_free(magnetic_space_group_type)
    type(c_ptr), value :: magnetic_space_group_type
end subroutine moyo_magnetic_space_group_type_free
```

Free a magnetic space-group type created by
`moyo_magnetic_space_group_type_new`.

### type `MoyoHallSymbolEntry`

- `hall_number` (`integer(c_int32_t)`) -- Hall number.
- `number` (`integer(c_int32_t)`) -- ITA space-group number.
- `arithmetic_number` (`integer(c_int32_t)`) -- arithmetic crystal class
  number.
- `setting` (`type(c_ptr)`) -- setting label C string; use `moyo_to_string`.
- `hall_symbol` (`type(c_ptr)`) -- Hall symbol C string.
- `hm_short` (`type(c_ptr)`) -- short Hermann-Mauguin symbol C string.
- `hm_full` (`type(c_ptr)`) -- full Hermann-Mauguin symbol C string.
- `centering` (`type(c_ptr)`) -- centering type C string.

### type `MoyoLayerHallSymbolEntry`

The components mirror `MoyoHallSymbolEntry`:

- `hall_number` (`integer(c_int32_t)`) -- layer Hall number.
- `number` (`integer(c_int32_t)`) -- layer-group number.
- `arithmetic_number` (`integer(c_int32_t)`) -- arithmetic crystal class
  number.
- `setting` (`type(c_ptr)`) -- setting label C string.
- `hall_symbol` (`type(c_ptr)`) -- Hall symbol C string.
- `hm_short` (`type(c_ptr)`) -- short Hermann-Mauguin symbol C string.
- `hm_full` (`type(c_ptr)`) -- full Hermann-Mauguin symbol C string.
- `centering` (`type(c_ptr)`) -- centering type C string.

### type `MoyoSpaceGroupType`

- `number` (`integer(c_int32_t)`) -- ITA space-group number.
- `hm_short` (`type(c_ptr)`) -- short Hermann-Mauguin symbol C string.
- `hm_full` (`type(c_ptr)`) -- full Hermann-Mauguin symbol C string.
- `arithmetic_number` (`integer(c_int32_t)`) -- arithmetic crystal class
  number.
- `arithmetic_symbol` (`type(c_ptr)`) -- arithmetic crystal class symbol C
  string.
- `geometric_crystal_class` (`type(c_ptr)`) -- geometric crystal class C
  string.
- `crystal_system` (`type(c_ptr)`) -- crystal system C string.
- `bravais_class` (`type(c_ptr)`) -- Bravais class C string.
- `lattice_system` (`type(c_ptr)`) -- lattice system C string.
- `crystal_family` (`type(c_ptr)`) -- crystal family C string.

### type `MoyoLayerGroupType`

- `number` (`integer(c_int32_t)`) -- layer-group number.
- `hm_short` (`type(c_ptr)`) -- short Hermann-Mauguin symbol C string.
- `hm_full` (`type(c_ptr)`) -- full Hermann-Mauguin symbol C string.
- `arithmetic_number` (`integer(c_int32_t)`) -- arithmetic crystal class
  number.
- `arithmetic_symbol` (`type(c_ptr)`) -- arithmetic crystal class symbol C
  string.
- `geometric_crystal_class` (`type(c_ptr)`) -- geometric crystal class C
  string.
- `bravais_class` (`type(c_ptr)`) -- Bravais class C string.
- `lattice_system` (`type(c_ptr)`) -- lattice system C string.

### type `MoyoMagneticSpaceGroupType`

- `uni_number` (`integer(c_int32_t)`) -- UNI number.
- `litvin_number` (`integer(c_int32_t)`) -- Litvin number.
- `bns_number` (`type(c_ptr)`) -- BNS number C string; use `moyo_to_string`.
- `og_number` (`type(c_ptr)`) -- OG number C string.
- `number` (`integer(c_int32_t)`) -- ITA space-group number of the family.
- `construct_type` (`integer(c_int32_t)`) -- magnetic space-group construction
  type.

## Group identification from primitive operations

```fortran
function moyo_point_group_new(prim_rotations, num_operations, basis) &
    result(point_group)
    type(c_ptr), value :: prim_rotations
    integer(c_int32_t), value :: num_operations
    type(c_ptr), value :: basis
    type(c_ptr) :: point_group
end function moyo_point_group_new
```

Identify the point group of `num_operations` primitive rotations.
`prim_rotations` points to the rotation matrices (for example `ops%rotations`,
or `c_loc` of an `integer(c_int32_t)` array whose `(:, :, k)` slices are the
transposed matrices). Pass `c_null_ptr` as `basis` to assume an identity basis,
or `c_loc` of a `real(c_double) :: basis(3, 3)` whose columns are the basis
vectors. Returns a pointer to `MoyoPointGroup`, or NULL on failure. Free with
`moyo_point_group_free`.

```fortran
subroutine moyo_point_group_free(point_group)
    type(c_ptr), value :: point_group
end subroutine moyo_point_group_free
```

Free a point group created by `moyo_point_group_new`.

```fortran
function moyo_space_group_new(prim_rotations, prim_translations, num_operations, &
                              basis, setting, hall_number, epsilon) &
    result(space_group)
    type(c_ptr), value :: prim_rotations
    type(c_ptr), value :: prim_translations
    integer(c_int32_t), value :: num_operations
    type(c_ptr), value :: basis
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    real(c_double), value :: epsilon
    type(c_ptr) :: space_group
end function moyo_space_group_new
```

Identify the space group from `num_operations` primitive rotations and
translations. Pass `c_null_ptr` as `basis` to assume an identity basis and a
negative `epsilon` to use the default tolerance (1e-4). `setting` is one of the
`MOYO_SETTING_*` constants and `hall_number` is only used with
`MOYO_SETTING_HALL_NUMBER`. Returns a pointer to `MoyoSpaceGroup`, or NULL on
failure. Free with `moyo_space_group_free`.

```fortran
subroutine moyo_space_group_free(space_group)
    type(c_ptr), value :: space_group
end subroutine moyo_space_group_free
```

Free a space group created by `moyo_space_group_new`.

```fortran
function moyo_layer_group_new(prim_rotations, prim_translations, num_operations, &
                              basis, setting, hall_number, epsilon) &
    result(layer_group)
    type(c_ptr), value :: prim_rotations
    type(c_ptr), value :: prim_translations
    integer(c_int32_t), value :: num_operations
    type(c_ptr), value :: basis
    integer(c_int32_t), value :: setting
    integer(c_int32_t), value :: hall_number
    real(c_double), value :: epsilon
    type(c_ptr) :: layer_group
end function moyo_layer_group_new
```

Identify the layer group from `num_operations` primitive layer-cell rotations
and translations. Pass `c_null_ptr` as `basis` to assume an identity basis and
a negative `epsilon` to use the default tolerance (1e-4). `setting` is one of
the `MOYO_LAYER_SETTING_*` constants and `hall_number` is only used with
`MOYO_LAYER_SETTING_HALL_NUMBER`. Returns a pointer to `MoyoLayerGroup`, or NULL
on failure. Free with `moyo_layer_group_free`.

```fortran
subroutine moyo_layer_group_free(layer_group)
    type(c_ptr), value :: layer_group
end subroutine moyo_layer_group_free
```

Free a layer group created by `moyo_layer_group_new`.

```fortran
function moyo_magnetic_space_group_new(prim_rotations, prim_translations, &
                                       prim_time_reversals, num_operations, &
                                       basis, epsilon) &
    result(magnetic_space_group)
    type(c_ptr), value :: prim_rotations
    type(c_ptr), value :: prim_translations
    type(c_ptr), value :: prim_time_reversals
    integer(c_int32_t), value :: num_operations
    type(c_ptr), value :: basis
    real(c_double), value :: epsilon
    type(c_ptr) :: magnetic_space_group
end function moyo_magnetic_space_group_new
```

Identify the magnetic space group from `num_operations` primitive magnetic
operations. Pass `c_null_ptr` as `basis` to assume an identity basis and a
negative `epsilon` to use the default tolerance (1e-4). Returns a pointer to
`MoyoMagneticSpaceGroup`, or NULL on failure. Free with
`moyo_magnetic_space_group_free`.

```fortran
subroutine moyo_magnetic_space_group_free(magnetic_space_group)
    type(c_ptr), value :: magnetic_space_group
end subroutine moyo_magnetic_space_group_free
```

Free a magnetic space group created by `moyo_magnetic_space_group_new`.

### type `MoyoPointGroup`

- `arithmetic_number` (`integer(c_int32_t)`) -- arithmetic crystal class
  number.
- `prim_trans_mat` (`integer(c_int32_t)(3, 3)`) -- transformation matrix
  (transpose view of the row-major C matrix).

### type `MoyoSpaceGroup`

- `number` (`integer(c_int32_t)`) -- ITA space-group number.
- `hall_number` (`integer(c_int32_t)`) -- Hall number.
- `linear` (`integer(c_int32_t)(3, 3)`) -- linear transformation (transpose
  view of the row-major C matrix).
- `origin_shift` (`real(c_double)(3)`) -- origin shift.

### type `MoyoLayerGroup`

- `number` (`integer(c_int32_t)`) -- layer-group number.
- `hall_number` (`integer(c_int32_t)`) -- layer Hall number.
- `linear` (`integer(c_int32_t)(3, 3)`) -- linear transformation (transpose
  view).
- `origin_shift` (`real(c_double)(3)`) -- origin shift.

### type `MoyoMagneticSpaceGroup`

- `uni_number` (`integer(c_int32_t)`) -- UNI number.
- `linear` (`integer(c_int32_t)(3, 3)`) -- linear transformation (transpose
  view).
- `origin_shift` (`real(c_double)(3)`) -- origin shift.

## Shared cell types

### type `MoyoCell`

- `basis` (`real(c_double)(3, 3)`) -- lattice basis; `basis(:, i)` is the `i`th
  basis vector.
- `positions` (`type(c_ptr)`) -- points to `num_atoms` fractional coordinate
  triplets; map with shape `[3, int(cell%num_atoms)]`.
- `numbers` (`type(c_ptr)`) -- points to `num_atoms` int32 atomic numbers.
- `num_atoms` (`integer(c_int32_t)`) -- number of atoms.

### type `MoyoCollinearMagneticCell`

- `cell` (`type(MoyoCell)`) -- the underlying crystal structure.
- `magnetic_moments` (`type(c_ptr)`) -- points to `cell%num_atoms` scalar
  moments.

### type `MoyoNonCollinearMagneticCell`

- `cell` (`type(MoyoCell)`) -- the underlying crystal structure.
- `magnetic_moments` (`type(c_ptr)`) -- points to `cell%num_atoms` Cartesian
  moment vectors; map with shape `[3, int(cell%num_atoms)]`.

## String helper

```fortran
function moyo_to_string(c_str) result(f_str)
    type(c_ptr), intent(in) :: c_str
    character(len=:), allocatable :: f_str
end function moyo_to_string
```

Convert a NUL-terminated C string (`type(c_ptr)`) into a Fortran `character`
value. Returns an empty string for a NULL pointer.

## Example

A compact space-group dataset example for an hcp structure:

```fortran
program example
    use moyo
    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: dataset_ptr
    type(MoyoDataset), pointer :: dataset
    real(c_double) :: basis(3, 3)
    real(c_double) :: positions(3, 2)
    integer(c_int32_t) :: numbers(2)
    real(c_double) :: a, c

    a = 3.17d0
    c = 5.14d0
    basis(:, 1) = [a, 0.0d0, 0.0d0]
    basis(:, 2) = [-a / 2.0d0, a * sqrt(3.0d0) / 2.0d0, 0.0d0]
    basis(:, 3) = [0.0d0, 0.0d0, c]
    positions(:, 1) = [1.0d0 / 3.0d0, 2.0d0 / 3.0d0, 1.0d0 / 4.0d0]
    positions(:, 2) = [2.0d0 / 3.0d0, 1.0d0 / 3.0d0, 3.0d0 / 4.0d0]
    numbers = [0_c_int32_t, 0_c_int32_t]

    dataset_ptr = moyo_dataset_new(basis, positions, numbers, 2_c_int32_t, &
                                   1d-4, -1d0, MOYO_SETTING_SPGLIB, &
                                   -1_c_int32_t, .true._c_bool)
    if (.not. c_associated(dataset_ptr)) error stop "moyo_dataset_new failed"

    call c_f_pointer(dataset_ptr, dataset)
    print *, "number:    ", dataset%number
    print *, "hm_symbol: ", moyo_to_string(dataset%hm_symbol)

    call moyo_dataset_free(dataset_ptr)
end program example
```
