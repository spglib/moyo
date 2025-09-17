cmake_minimum_required(VERSION 3.22)
project(example_moyoc LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    moyoc
    URL      https://github.com/spglib/moyo/releases/download/v0.4.7/moyoc-0.4.7-aarch64-apple-darwin.tar.gz
    URL_HASH SHA256=c792df87fabf620ca5e473fb6db99d345656da94031426292173a589a089d72f
)
FetchContent_MakeAvailable(moyoc)

add_library(moyoc SHARED IMPORTED)
set_target_properties(moyoc PROPERTIES
    IMPORTED_LOCATION "${moyoc_SOURCE_DIR}/lib/libmoyoc.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "${moyoc_SOURCE_DIR}/include"
    BUILD_RPATH "${moyoc_SOURCE_DIR}/lib"
)

add_executable(example example.c)
target_link_libraries(example PRIVATE moyoc)
