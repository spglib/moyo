function(add_moyoc_test name)
    add_executable(${name} ${name}.c)
    target_link_libraries(${name} PRIVATE moyo::moyoc)
    if(MSVC)
        # Re-enable assert() in Release configurations
        target_compile_options(${name} PRIVATE /UNDEBUG)
    else()
        target_compile_options(${name} PRIVATE
            # Re-enable assert() in Release configurations
            -UNDEBUG
            -Wall -Wextra -fsanitize=address -fno-omit-frame-pointer)
        target_link_options(${name} PRIVATE -fsanitize=address)
        target_link_libraries(${name} PRIVATE m)
    endif()

    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_moyoc_test(test_moyoc_data_entries)
add_moyoc_test(test_moyoc_dataset)
add_moyoc_test(test_moyoc_identify)
add_moyoc_test(test_moyoc_layer_dataset)
add_moyoc_test(test_moyoc_magnetic_dataset)
add_moyoc_test(test_moyoc_operations)
