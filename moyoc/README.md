# moyoc

C bindings for [moyo](https://github.com/spglib/moyo), built with CMake and
[corrosion](https://github.com/corrosion-rs/corrosion).

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
See [examples](examples) for a complete project.

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

## API

See the generated `moyoc.h` for the full API. The main entry points are:

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
```
