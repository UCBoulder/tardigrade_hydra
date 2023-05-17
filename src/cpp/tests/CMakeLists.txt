set(TEST_NAME "test_${PROJECT_NAME}")
add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
target_compile_options(${TEST_NAME} PRIVATE)
target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME} ${PROJECT_LINK_LIBRARIES})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${TEST_NAME} PUBLIC
                               ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                               ${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${abaqus_tools_SOURCE_DIR}/${CPP_SRC_PATH})
endif()

foreach(support_module ${INTERNAL_SUPPORT_LIBRARIES})

    set(TEST_NAME "test_${support_module}")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    target_compile_options(${TEST_NAME} PRIVATE "-lrt")
    target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME})

    # Local builds of upstream projects require local include paths
    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(${TEST_NAME} PUBLIC
                                   ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                                   ${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${stress_tools_SOURCE_DIR}/${CPP_SRC_PATH})
    endif()

endforeach(support_module)
