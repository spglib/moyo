add_executable(test_moyoc_dataset test_moyoc_dataset.c)
target_link_libraries(test_moyoc_dataset PRIVATE moyo::moyoc)
target_compile_features(test_moyoc_dataset PRIVATE c_std_11)
if(NOT MSVC)
    target_compile_options(test_moyoc_dataset PRIVATE
        -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(test_moyoc_dataset PRIVATE -fsanitize=address)
    target_link_libraries(test_moyoc_dataset PRIVATE m)
endif()

# Locate libmoyoc at runtime when linking the shared library
if(BUILD_SHARED_LIBS AND TARGET moyoc-shared)
    set_property(TARGET test_moyoc_dataset
        PROPERTY BUILD_RPATH "$<TARGET_FILE_DIR:moyoc-shared>")
endif()

add_test(NAME test_moyoc_dataset COMMAND test_moyoc_dataset)
