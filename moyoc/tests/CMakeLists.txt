add_executable(test_moyoc_dataset test_moyoc_dataset.c)
target_link_libraries(test_moyoc_dataset PRIVATE moyo::moyoc)
if(NOT MSVC)
    target_compile_options(test_moyoc_dataset PRIVATE
        -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(test_moyoc_dataset PRIVATE -fsanitize=address)
    target_link_libraries(test_moyoc_dataset PRIVATE m)
endif()

add_test(NAME test_moyoc_dataset COMMAND test_moyoc_dataset)
