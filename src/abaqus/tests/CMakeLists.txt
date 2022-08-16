if(cmake_build_type_lower STREQUAL "conda-test")
    set(umat_file_string ${installed_umat_file})
else()
    set(umat_file_string $<TARGET_FILE:${UMAT}>)
endif()
# Add abaqus cpp_stub_umat integration tests that run abaqus against input files
foreach(INPUT_FILE ${ABAQUS_INPUT_FILES})
    add_test(NAME ${ABAQUS_TEST_NAME}[${INPUT_FILE}]
             COMMAND ${BASH_PROGRAM} ${ABAQUS_INTEGRATION_TEST} ${ABAQUS_PROGRAM} ../${INPUT_FILE} ${umat_file_string}
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}/tests
            )
endforeach()
