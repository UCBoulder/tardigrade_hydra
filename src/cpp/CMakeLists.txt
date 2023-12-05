add_library (${PROJECT_NAME} SHARED ${PROJECT_SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${PROJECT_NAME}.h)
target_link_libraries(${PROJECT_NAME} tardigrade_constitutive_tools tardigrade_stress_tools ${PROJECT_LINK_LIBRARIES})
target_compile_options(${PROJECT_NAME} PUBLIC)

# Abaqus UMAT interface
add_library(${UMAT} SHARED "${UMAT}.cpp" "${UMAT}.h")
set_target_properties(${UMAT} PROPERTIES
                              PUBLIC_HEADER ${UMAT}.h
                              PREFIX ""
                              SUFFIX ".o")
target_compile_options(${UMAT} PUBLIC)
target_link_libraries(${UMAT} PUBLIC ${PROJECT_NAME})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${PROJECT_NAME} PUBLIC
                               "${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${tardigrade_abaqus_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
    target_include_directories(${UMAT} PUBLIC
                               "${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
endif()

foreach(package ${PROJECT_NAME} ${UMAT})
    install(TARGETS ${package}
            EXPORT ${package}_Targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endforeach(package)
install(FILES ${PROJECT_PRIVATE_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
