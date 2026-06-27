# moyoc

C bindings for [moyo](https://github.com/spglib/moyo), built with CMake and
[corrosion](https://github.com/corrosion-rs/corrosion).

- [API reference](api/index.html)

## Requirements

- CMake >= 3.23
- Rust toolchain (cargo)
- [cbindgen](https://github.com/mozilla/cbindgen) (optional; built automatically if not installed)

## Using moyoc from a CMake project (recommended)

Build moyoc from source as part of your project with `FetchContent`:

```cmake
include(FetchContent)
FetchContent_Declare(
    moyo
    GIT_REPOSITORY https://github.com/spglib/moyo.git
    GIT_TAG v0.12.0 # Pin a release tag
    SOURCE_SUBDIR moyoc
)
FetchContent_MakeAvailable(moyo)

add_executable(example example.c)
target_link_libraries(example PRIVATE moyo::moyoc)
```

The `moyo::moyoc` target carries the generated `moyoc.h` include directory and links the
static library by default. Set `-DBUILD_SHARED_LIBS=ON` to link the shared library instead.
See [examples](https://github.com/spglib/moyo/tree/main/moyoc/examples) for a complete project.

## Standalone build

```shell
cmake -S . -B build
cmake --build build
ctest --test-dir build
cmake --install build --prefix <prefix>
```

This installs `include/moyoc.h`, `lib/libmoyoc.a`, and `lib/libmoyoc.{so,dylib}` for
non-CMake consumers. Prebuilt archives with this layout are attached to
[GitHub releases](https://github.com/spglib/moyo/releases).

From the repo root, `just c-build` and `just c-test` wrap the configure/build/test steps.

## API

See the [API reference](api/index.html) (generated with
Doxygen from `moyoc.h`; build it locally with `just c-docs`). The main entry
points are:

```c
// Space-group symmetry analysis
MoyoDataset *dataset = moyo_dataset_new(
    basis, positions, numbers, num_atoms,
    symprec, angle_tolerance, setting, hall_number, rotate_basis
);
// ... read dataset fields ...
moyo_dataset_free(dataset);

// Layer-group symmetry analysis
MoyoLayerDataset *layer_dataset = moyo_layer_dataset_new(
    basis, positions, numbers, num_atoms,
    symprec, angle_tolerance, layer_setting, hall_number, rotate_basis
);
// ... read layer_dataset fields ...
moyo_layer_dataset_free(layer_dataset);

// Magnetic space-group symmetry analysis (collinear moments)
MoyoCollinearMagneticDataset *mag_dataset = moyo_collinear_magnetic_dataset_new(
    basis, positions, numbers, magnetic_moments, num_atoms,
    symprec, angle_tolerance, mag_symprec, is_axial, rotate_basis
);
// ... read mag_dataset fields ...
moyo_collinear_magnetic_dataset_free(mag_dataset);

// Magnetic space-group symmetry analysis (non-collinear moments)
MoyoNonCollinearMagneticDataset *nc_mag_dataset = moyo_noncollinear_magnetic_dataset_new(
    basis, positions, numbers, cartesian_magnetic_moments, num_atoms,
    symprec, angle_tolerance, mag_symprec, is_axial, rotate_basis
);
// ... read nc_mag_dataset fields ...
moyo_noncollinear_magnetic_dataset_free(nc_mag_dataset);

// Symmetry operations from the database
MoyoOperations *operations = moyo_operations_from_number(
    number, setting, hall_number, primitive
);
moyo_operations_free(operations);

MoyoOperations *layer_operations = moyo_operations_from_layer_number(
    layer_number, layer_setting, layer_hall_number, primitive
);
moyo_operations_free(layer_operations);

MoyoMagneticOperations *magnetic_operations =
    moyo_magnetic_operations_from_uni_number(uni_number, primitive);
moyo_magnetic_operations_free(magnetic_operations);

// Group identification from primitive operations
MoyoPointGroup *point_group = moyo_point_group_new(
    prim_rotations, num_operations, basis /* or NULL */
);
moyo_point_group_free(point_group);

MoyoSpaceGroup *space_group = moyo_space_group_new(
    prim_rotations, prim_translations, num_operations,
    basis /* or NULL */, setting, hall_number, epsilon
);
moyo_space_group_free(space_group);

MoyoLayerGroup *layer_group = moyo_layer_group_new(
    prim_rotations, prim_translations, num_operations,
    basis /* or NULL */, layer_setting, layer_hall_number, epsilon
);
moyo_layer_group_free(layer_group);

MoyoMagneticSpaceGroup *magnetic_space_group = moyo_magnetic_space_group_new(
    prim_rotations, prim_translations, prim_time_reversals, num_operations,
    basis /* or NULL */, epsilon
);
moyo_magnetic_space_group_free(magnetic_space_group);

// Group-type metadata lookups
MoyoHallSymbolEntry *entry = moyo_hall_symbol_entry_new(hall_number);
moyo_hall_symbol_entry_free(entry);

MoyoLayerHallSymbolEntry *layer_entry = moyo_layer_hall_symbol_entry_new(layer_hall_number);
moyo_layer_hall_symbol_entry_free(layer_entry);

MoyoSpaceGroupType *space_group_type = moyo_space_group_type_new(number);
moyo_space_group_type_free(space_group_type);

MoyoLayerGroupType *layer_group_type = moyo_layer_group_type_new(layer_number);
moyo_layer_group_type_free(layer_group_type);

MoyoMagneticSpaceGroupType *magnetic_space_group_type =
    moyo_magnetic_space_group_type_new(uni_number);
moyo_magnetic_space_group_type_free(magnetic_space_group_type);
```

## Fortran interface

Configure with `-DMOYOC_BUILD_FORTRAN=ON` to build the Fortran interface and expose
the `moyo::moyof` CMake target (requires a Fortran compiler, e.g. gfortran). Link it
and `use moyo` from your Fortran code.

The interface mirrors the C API with a few conventions:

- Arrays are passed column-major: `basis(:, i)` is the i-th basis vector and
  `positions(:, i)` is the i-th site, i.e. a transpose view of the C row-major arrays.
- Constructors return a `type(c_ptr)` which is NULL on failure. Map it to the
  matching `bind(c)` derived type with `c_f_pointer`, then free it with the
  corresponding `*_free` routine.
- C strings are converted to Fortran with `moyo_to_string`.

```fortran
use moyo
use iso_c_binding, only: c_ptr, c_f_pointer, c_associated
type(c_ptr) :: handle
type(MoyoDataset), pointer :: dataset

handle = moyo_dataset_new(basis, positions, numbers, num_atoms, &
    symprec, angle_tolerance, setting, hall_number, rotate_basis)
if (.not. c_associated(handle)) error stop "moyo_dataset_new failed"
call c_f_pointer(handle, dataset)
! ... read dataset fields ...
call moyo_dataset_free(handle)
```

## How to cite moyoc

See the citation information in [the root README](https://github.com/spglib/moyo/blob/main/README.md)
