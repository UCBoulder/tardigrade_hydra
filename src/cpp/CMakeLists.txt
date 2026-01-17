add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCE_FILES})
set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES PUBLIC_HEADER ${PROJECT_NAME}.h
)
target_link_libraries(
    ${PROJECT_NAME}
    tardigrade_constitutive_tools
    tardigrade_stress_tools
    tardigrade_micromorphic_tools
    ${PROJECT_LINK_LIBRARIES}
)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${PROJECT_NAME}
        PUBLIC ${LOCAL_BUILD_INCLUDE_DIRECTORIES}
    )
endif()

foreach(package ${PROJECT_NAME})
    install(
        TARGETS ${package}
        EXPORT ${package}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    get_target_property(target_sources ${package} SOURCES)
    foreach(source ${target_sources})
        set_property(
            GLOBAL
            APPEND
            PROPERTY CLANG_FORMAT_SOURCE_FILES "${PROJECT_SOURCE_DIR}/${source}"
        )
    endforeach(source)
endforeach(package)

foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
    install(
        TARGETS ${package}
        EXPORT ${package}_Targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
    install(
        FILES ${package}.h ${package}.cpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    get_target_property(target_sources ${package} SOURCES)
    foreach(source ${target_sources})
        set_property(
            GLOBAL
            APPEND
            PROPERTY CLANG_FORMAT_SOURCE_FILES "${PROJECT_SOURCE_DIR}/${source}"
        )
    endforeach(source)
endforeach(package)

install(
    FILES ${PROJECT_PRIVATE_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

get_target_property(target_sources ${PROJECT_NAME} SOURCES)
foreach(source ${target_sources})
    set_property(
        GLOBAL
        APPEND
        PROPERTY
            CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}"
    )
endforeach(source)
