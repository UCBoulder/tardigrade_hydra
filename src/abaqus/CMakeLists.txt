# Explicitly set list of input files for integration testing
set(ABAQUS_INPUT_FILES
    "single_element_c3d8.inp"
   )
# Create variable for hardcoded integration shell script name
set(ABAQUS_INTEGRATION_TEST "test_abaqus.sh")
# Create integration test target and associated dependencies
set(ABAQUS_TEST_NAME "test_abaqus_integration")
add_custom_target(${ABAQUS_TEST_NAME} ALL
                  DEPENDS ${UMAT})
# Copy abaqus source directory to build directory
add_custom_command(TARGET ${ABAQUS_TEST_NAME} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       "${PROJECT_SOURCE_DIR}/${ABAQUS_SRC_PATH}" "${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}"
                   DEPENDS ${ABAQUS_INPUT_FILES} ${ABAQUS_INTEGRATION_TEST}
                   COMMENT "Preparing the ${PROJECT_NAME} Abaqus integration test files"
                  )

# Abaqus integration tests
if(BASH_PROGRAM AND ABAQUS_PROGRAM)
    add_subdirectory(./tests)
endif()
