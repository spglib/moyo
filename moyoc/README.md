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

From the repo root, `just c-build` and `just c-test` wrap the configure/build/test steps.

## API

See the [C API reference](docs/api.md) for the complete API (datasets,
database lookups, group-type metadata, and group identification), or the
generated `moyoc.h`. The main entry point is the space-group dataset:

```c
MoyoDataset *dataset = moyo_dataset_new(
    basis, positions, numbers, num_atoms,
    symprec, angle_tolerance, setting, hall_number, rotate_basis
);
// ... read dataset fields ...
moyo_dataset_free(dataset);
```

## Fortran interface

Configure with `-DMOYOC_BUILD_FORTRAN=ON` to build the Fortran interface and expose
the `moyo::moyof` CMake target (requires a Fortran compiler, e.g. gfortran). Link it
and `use moyo` from your Fortran code. See the
[Fortran API reference](docs/api_fortran.md) for the complete API.

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

See the citation information in [the root README](../README.md)
