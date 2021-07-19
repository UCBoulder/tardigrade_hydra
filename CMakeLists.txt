#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.14)
project(cpp_stub)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(ABAQUS_SRC_PATH "src/abaqus")
set(ABAQUS_TEST_PATH "${ABAQUS_SRC_PATH}/tests")
set(DOCUMENTATION_SRC_PATH "docs")
# Set some Abaqus integration variables
set(UMAT "cpp_stub_umat")
set(MISSING_INTEGRATION_DEPS_MESSAGE "Abaqus integration tests will be skipped.")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
# FIXME: Figure out why include() doesn't pick up CMAKE_MODULE_PATH correctly
include(${CMAKE_SRC_PATH}/GetVersionFromGitTag.cmake)
project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb")

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find bash (Required for abaqus integration tests)
find_program(BASH_PROGRAM bash)
if(BASH_PROGRAM)
    message(STATUS "Found bash: ${BASH_PROGRAM}")
else()
    message(STATUS "bash not found. ${MISSING_INTEGRATION_DEPS_MESSAGE}")
endif()

# Find Abaqus (Required for abaqus integration tests)
find_program(ABAQUS_PROGRAM
             NAMES abaqus abq2020 abq2019hf4 abq2019 abq2018hf11 abq2018 abq2017hf2 abq2017
             HINTS /apps/abaqus/Commands
            )
if(ABAQUS_PROGRAM)
    message(STATUS "Found Abaqus: ${ABAQUS_PROGRAM}")
else()
    message(STATUS "Abaqus not found. ${MISSING_INTEGRATION_DEPS_MESSAGE}")
endif()

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
endif()

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "error_tools" "vector_tools" "abaqus_tools")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    find_package(${package} ${upstream_required} CONFIG)
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        # TODO: Decide how to handle master/dev build differences. bash env variables for GIT TAG?
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY ssh://git@re-git.lanl.gov:10022/aea/material-models/${package}.git
            GIT_TAG        origin/dev
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
include_directories(${CPP_SRC_PATH})
add_subdirectory(${CPP_SRC_PATH})

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add c++ tests and docs
    add_subdirectory(${CPP_TEST_PATH})
    add_subdirectory(${DOCUMENTATION_SRC_PATH})
    add_subdirectory(${ABAQUS_SRC_PATH})
    # Abaqus integration tests
    if(BASH_PROGRAM AND ABAQUS_PROGRAM)
        add_subdirectory(${ABAQUS_TEST_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
foreach(package ${PROJECT_NAME} ${UMAT})
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${package}ConfigVersion.cmake"
                                     VERSION ${PROJECT_VERSION}
                                     COMPATIBILITY SameMajorVersion)
    configure_package_config_file(
      "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
      "${PROJECT_BINARY_DIR}/${package}Config.cmake"
      INSTALL_DESTINATION
      ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

    # CMake won't build the targets for local builds of upstream projects
    if(cmake_build_type_lower STREQUAL release)
        install(EXPORT ${package}_Targets
                FILE ${package}Targets.cmake
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
    endif()

    install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                  "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
endforeach(package)
