cmake_minimum_required(VERSION 3.23)

# Read the version from the workspace manifest
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../Cargo.toml" _moyo_workspace_manifest)
string(REGEX MATCH "\nversion = \"([0-9]+\\.[0-9]+\\.[0-9]+)\"" _ "${_moyo_workspace_manifest}")
if(NOT CMAKE_MATCH_1)
    message(FATAL_ERROR "Failed to parse the version from the workspace Cargo.toml")
endif()

project(moyoc VERSION "${CMAKE_MATCH_1}" LANGUAGES C)

option(MOYOC_BUILD_TESTING "Build moyoc C tests" ${PROJECT_IS_TOP_LEVEL})
option(MOYOC_INSTALL "Generate install rules for moyoc" ${PROJECT_IS_TOP_LEVEL})

if(PROJECT_IS_TOP_LEVEL)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    endif()
    # Corrosion copies the cargo-built libraries into these directories
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

find_package(Corrosion 0.6 QUIET)
if(NOT Corrosion_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        Corrosion
        GIT_REPOSITORY https://github.com/corrosion-rs/corrosion.git
        GIT_TAG v0.6.1
    )
    FetchContent_MakeAvailable(Corrosion)
endif()

# Creates an INTERFACE target `moyoc` linking `moyoc-static` or `moyoc-shared`
# depending on BUILD_SHARED_LIBS
corrosion_import_crate(
    MANIFEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml"
    CRATES moyoc
)

# Give the shared library a proper install name / soname so that installed
# copies are usable
if(Rust_CARGO_TARGET_OS STREQUAL "darwin")
    corrosion_add_target_local_rustflags(moyoc "-Clink-arg=-Wl,-install_name,@rpath/libmoyoc.dylib")
elseif(Rust_CARGO_TARGET_OS STREQUAL "linux")
    corrosion_add_target_local_rustflags(moyoc "-Clink-arg=-Wl,-soname,libmoyoc.so")
endif()

# Generate moyoc.h with cbindgen and attach it to the `moyoc` target
corrosion_experimental_cbindgen(
    TARGET moyoc
    HEADER_NAME moyoc.h
)
set(MOYOC_GENERATED_HEADER
    "${CMAKE_CURRENT_BINARY_DIR}/corrosion_generated/cbindgen/moyoc/include/moyoc.h")

add_library(moyo::moyoc ALIAS moyoc)

if(MOYOC_BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

if(MOYOC_INSTALL)
    include(GNUInstallDirs)
    install(FILES "${MOYOC_GENERATED_HEADER}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(FILES "$<TARGET_FILE:moyoc-static>" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    if(Rust_CARGO_TARGET_OS STREQUAL "windows")
        install(FILES "$<TARGET_FILE:moyoc-shared>" DESTINATION "${CMAKE_INSTALL_BINDIR}")
        install(FILES "$<TARGET_LINKER_FILE:moyoc-shared>" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    else()
        install(FILES "$<TARGET_FILE:moyoc-shared>" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    endif()
endif()
