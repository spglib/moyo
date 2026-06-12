cmake_minimum_required(VERSION 3.23)
project(example_moyoc LANGUAGES C)

include(FetchContent)
FetchContent_Declare(
    moyo
    GIT_REPOSITORY https://github.com/spglib/moyo.git
    GIT_TAG main # Pin a release tag (e.g. v0.12.0) in real projects
    SOURCE_SUBDIR moyoc
)
FetchContent_MakeAvailable(moyo)

add_executable(example example.c)
target_link_libraries(example PRIVATE moyo::moyoc)
if(NOT MSVC)
    target_link_libraries(example PRIVATE m)
endif()
