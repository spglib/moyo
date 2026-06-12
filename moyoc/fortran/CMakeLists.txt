add_library(moyof STATIC moyof.f90)
target_link_libraries(moyof PUBLIC moyo::moyoc)

# Expose the generated `moyo.mod` to consumers of `moyof`
set_target_properties(moyof PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
target_include_directories(moyof INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

add_library(moyo::moyof ALIAS moyof)

if(MOYOC_BUILD_TESTING)
    add_executable(test_moyof_dataset tests/test_moyof_dataset.f90)
    target_link_libraries(test_moyof_dataset PRIVATE moyo::moyof)
    add_test(NAME test_moyof_dataset COMMAND test_moyof_dataset)
endif()

if(MOYOC_INSTALL)
    include(GNUInstallDirs)
    install(FILES "$<TARGET_FILE:moyof>" DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    # .mod files are compiler-specific; they are only usable by the same
    # Fortran compiler (and version) that produced them.
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/moyo.mod"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif()
