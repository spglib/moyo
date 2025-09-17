cmake_minimum_required(VERSION 3.22)
project(example_moyoc LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    moyoc
    URL      https://github.com/spglib/moyo/releases/download/v0.4.8/moyoc-0.4.8-aarch64-apple-darwin.tar.gz
    URL_HASH SHA256=2eccd32ec7f4c33088b2d0d2c18731d3b920a1a57948948a4b97b0bdabedbd8f
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
