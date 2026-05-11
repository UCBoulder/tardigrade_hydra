#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.30)
cmake_policy(SET CMP0110 NEW)
project(tardigrade_hydra LANGUAGES CXX)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")

# Add a flag for whether the python bindings should be built or not
set(TARDIGRADE_HYDRA_BUILD_PYTHON_BINDINGS
    ON
    CACHE BOOL
    "Flag for whether the python bindings should be built for hydra"
)
set(TARDIGRADE_HYDRA_USE_LIBXSMM OFF CACHE BOOL "Flag for whether to use libxsmm for matrix math")

# Add a flag for if a full build of all tardigrade repositories should be performed
set(TARDIGRADE_FULL_BUILD
    OFF
    CACHE BOOL
    "Flag for whether a full build of Tardigrade should be performed (i.e., all repos pulled from git and built)"
)

if(${TARDIGRADE_ERROR_TOOLS_OPT})
    add_compile_definitions(TARDIGRADE_ERROR_TOOLS_OPT)
    message(WARNING "BUILDING OPTIMIZED ERROR TOOLS. NO ERRORS WILL BE CAUGHT")
endif()

# Set the internal support libraries
set(INTERNAL_SUPPORT_LIBRARIES
    "tardigrade_hydraMicromorphic"
    "tardigrade_hydraLinearElasticity"
    "tardigrade_hydraLinearTestMaterial"
    "tardigrade_hydraLinearViscoelasticity"
    "tardigrade_hydraThermalExpansion"
    "tardigrade_hydraPerzynaViscoplasticity"
    "tardigrade_hydraMicromorphicLinearElasticity"
    "tardigrade_hydraPerzynaViscodamage"
    "tardigrade_hydraMicromorphicDruckerPragerPlasticity"
    "tardigrade_hydraPerzynaJ2Viscoplasticity"
    "tardigrade_hydraMassChange"
    "tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization"
    "tardigrade_hydraFourierHeatConduction"
    "tardigrade_hydraLinearInternalEnergy"
    "tardigrade_hydraMicromorphicRadialReturnLinearElasticity"
    "tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity"
    "tardigrade_hydraDOFVelocityGradientDeformation"
    "tardigrade_SolverBase"
    "tardigrade_IterativeSolverBase"
    "tardigrade_RelaxedSolverBase"
    "tardigrade_RelaxedSolver"
    "tardigrade_SubcyclerSolver"
    "tardigrade_NewtonSolver"
    "tardigrade_StepDampingBase"
    "tardigrade_GradientDamping"
    "tardigrade_ArmijoGradientDamping"
    "tardigrade_SolverStepBase"
    "tardigrade_SolverStep"
    "tardigrade_ArmijoLineSearchStep"
    "tardigrade_PreconditionerBase"
    "tardigrade_MaxRowPreconditioner"
    "tardigrade_TrialStepBase"
    "tardigrade_NonlinearStepBase"
    "tardigrade_NewtonStep"
    "tardigrade_SequentialQuadraticProgrammingStep"
    "tardigrade_LevenbergMarquardtStep"
    "tardigrade_DeformationBase"
    "tardigrade_HyperelasticBase"
    "tardigrade_NeoHookianStrainEnergy"
    "tardigrade_MooneyRivlinStrainEnergy"
    "tardigrade_CHIPFoamStrainEnergy"
)
set(ADDITIONAL_HEADER_ONLY_LIBRARIES
    "tardigrade_DeformationDecompositionBase"
    "tardigrade_SetDataStorage"
    "tardigrade_MatrixMap"
    "tardigrade_CoreDefinitions"
    "tardigrade_CustomErrors"
    "tardigrade_ResidualBase"
    "tardigrade_DeformationEvolutionBase"
)
set(PROJECT_SOURCE_FILES ${PROJECT_NAME}.cpp ${PROJECT_NAME}.h ${PROJECT_NAME}.tpp)
set(PROJECT_PRIVATE_HEADERS "")
foreach(package ${INTERNAL_SUPPORT_LIBRARIES})
    set(PROJECT_SOURCE_FILES ${PROJECT_SOURCE_FILES} ${package}.cpp ${package}.h ${package}.tpp)
    set(PROJECT_PRIVATE_HEADERS ${PROJECT_PRIVATE_HEADERS} ${package}.h)
endforeach(package)
if(TARDIGRADE_HYDRA_USE_LIBXSMM)
    message(STATUS "Using libxsmm for matrix math")
    find_library(libxsmm_DIR xsmm REQUIRED PATHS "$ENV{CONDA_PREFIX}/lib")
    add_definitions(-DTARDIGRADE_HYDRA_USE_LLXSMM)
    set(PROJECT_LINK_LIBRARIES ${PROJECT_LINK_LIBRARIES} "xsmm" "blas")
endif()
# Upstream libraries for each sub-library
set(TARDIGRADE_HYDRAMICROMORPHIC_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools")
set(TARDIGRADE_HYDRALINEARELASTICITY_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools")
set(TARDIGRADE_HYDRALINEARVISCOELASTICITY_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools" "tardigrade_stress_tools")
set(TARDIGRADE_HYDRATHERMALEXPANSION_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools")
set(TARDIGRADE_HYDRAPERYZNAVISCOPLASTICITY_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools" "tardigrade_stress_tools")
set(TARDIGRADE_HYDRAPERYZNAVISCODAMAGE_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools" "tardigrade_stress_tools")
set(TARDIGRADE_HYDRAPERYZNAJ2VISCOPLASTICITY_UPSTREAM_LIBRARIES
    "tardigrade_constitutive_tools"
    "tardigrade_stress_tools"
)
set(TARDIGRADE_HYDRAMASSCHANGE_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools" "tardigrade_stress_tools")
set(TARDIGRADE_HYDRAMICROMORPHICLINEARELASTICITY_UPSTREAM_LIBRARIES
    "tardigrade_constitutive_tools"
    "tardigrade_micromorphic_tools"
)
set(TARDIGRADE_HYDRAMICROMORPHICLINEARTESTMATERIAL_UPSTREAM_LIBRARIES
    "tardigrade_constitutive_tools"
    "tardigrade_micromorphic_tools"
)
set(TARDIGRADE_HYDRAMICROMORPHICDRUCKERPRAGERPLASTICITY_UPSTREAM_LIBRARIES
    "tardigrade_constitutive_tools"
    "tardigrade_micromorphic_tools"
)
set(TARDIGRADE_HYDRAMASSCHANGEDOF_UPSTREAM_LIBRARIES "tardigrade_constitutive_tools")

# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
endif()

if(${not_conda_test} STREQUAL "true")
    foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
        add_library(
            ${package}
            INTERFACE
            "${CPP_SRC_PATH}/${package}.cpp"
            "${CPP_SRC_PATH}/${package}.h"
            "${CPP_SRC_PATH}/${package}.tpp"
        )
    endforeach(package)
endif()

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()

# Get version number from Git
set(VERSION_UPDATE_FROM_GIT True)
if(${not_conda_test} STREQUAL "true")
    if(DEFINED ENV{PREFIX})
        set(Python_ROOT_DIR "$ENV{BUILD_PREFIX}/bin")
    endif()
    set(Python_FIND_STRATEGY LOCATION)
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(
        COMMAND ${Python_EXECUTABLE} -m setuptools_scm
        OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(${${PROJECT_NAME}_VERSION_STRING_FULL} STREQUAL "")
        set(${PROJECT_NAME}_VERSION 0.0.0)
    else()
        string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    endif()
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===

#=============================================================================================== FIND PROGRAMS ===

find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format REQUIRED)
set_property(GLOBAL PROPERTY CLANG_FORMAT_SOURCE_FILES)

#=============================================================================================== ADD PROJECT TARGETS ===

# Find eigen
find_package(Eigen3 3.3 NO_MODULE)
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3 (less than version 3.4.1): ${EIGEN3_INCLUDE_DIR}")
else()
    find_package(Eigen3 3.4...5 NO_MODULE)
    if(EIGEN3_FOUND)
        message(STATUS "Found Eigen3 (less than version 5): ${EIGEN3_INCLUDE_DIR}")
    else()
        find_package(Eigen3 5 REQUIRED NO_MODULE)
        message(STATUS "Found Eigen3 (version 5): ${EIGEN3_INCLUDE_DIR}")
    endif()
endif()

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages
    "tardigrade_error_tools"
    "tardigrade_vector_tools"
    "tardigrade_constitutive_tools"
    "tardigrade_stress_tools"
    "tardigrade_micromorphic_tools"
)
if(${not_conda_test} STREQUAL "false")
    # During conda-build testing, we must find the installed project files as if they were an external project
    set(upstream_packages ${upstream_packages} ${PROJECT_NAME})
endif()
set(LOCAL_BUILD_INCLUDE_DIRECTORIES "")
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS
        ${TARDIGRADE_HYDRA_BUILD_PYTHON_BINDINGS}
        CACHE INTERNAL
        "Setting ${package}'s python binding flag to the global value"
    )
    if(TARDIGRADE_FULL_BUILD)
        message("Not attempting to find ${package}")
    else()
        find_package(${package} ${upstream_required} CONFIG)
    endif()
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        if(NOT DEFINED ${package_upper}_FETCHCONTENT_VERSION)
            set(${package_upper}_FETCHCONTENT_VERSION "origin/dev")
        endif()
        message("${package_upper} is being built with version ${${package_upper}_FETCHCONTENT_VERSION}")
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY https://github.com/UCBoulder/${package}.git
            GIT_TAG ${${package_upper}_FETCHCONTENT_VERSION}
        )
        FetchContent_MakeAvailable(${package})
        set(LOCAL_BUILD_INCLUDE_DIRECTORIES ${LOCAL_BUILD_INCLUDE_DIRECTORIES} ${${package}_SOURCE_DIR}/${CPP_SRC_PATH})
    endif()
endforeach(package)

foreach(package ${upstream_packages})
    if(NOT cmake_build_type_lower STREQUAL "release")
        string(TOLOWER "${package}" lower_package)
        foreach(sub_package ${${upper_package}_UPSTREAM_LIBRARIES})
            target_include_directories(${package} PUBLIC ${${sub_package}_SOURCE_DIR}/${CPP_SRC_PATH})
        endforeach(sub_package)
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
if(${not_conda_test} STREQUAL "true")
    include_directories(${CPP_SRC_PATH})
    add_subdirectory(${CPP_SRC_PATH})
endif()

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add c++ tests and docs
    add_subdirectory(${CPP_TEST_PATH})
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ===
if(${not_conda_test} STREQUAL "true")
    foreach(package ${PROJECT_NAME})
        include(CMakePackageConfigHelpers)
        write_basic_package_version_file(
            "${package}ConfigVersion.cmake"
            VERSION ${PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
        )
        configure_package_config_file(
            "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
            "${PROJECT_BINARY_DIR}/${package}Config.cmake"
            INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
        )

        # CMake won't build the targets for local builds of upstream projects
        if(cmake_build_type_lower STREQUAL release)
            install(
                EXPORT ${package}_Targets
                FILE ${package}Targets.cmake
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
            )
        endif()

        install(
            FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake" "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake
        )
    endforeach(package)
    install(
        FILES ${CPP_SRC_PATH}/${PROJECT_NAME}.cpp ${CPP_SRC_PATH}/${PROJECT_NAME}.tpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
    foreach(package ${INTERNAL_SUPPORT_LIBRARIES})
        install(FILES ${CPP_SRC_PATH}/${package}.cpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
        install(FILES ${CPP_SRC_PATH}/${package}.tpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endforeach(package)
endif()

#================================================================================== SETUP CLANG FORMATTING TARGETS ===

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    get_property(CLANG_FORMAT_SOURCE_FILES_VALUE GLOBAL PROPERTY CLANG_FORMAT_SOURCE_FILES)

    add_custom_target(
        cpp-format-check
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -style=file --dry-run -Werror -i ${CLANG_FORMAT_SOURCE_FILES_VALUE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Check the C++ code format with clang-format"
    )
    set_target_properties(cpp-format-check PROPERTIES EXCLUDE_FROM_ALL TRUE)

    add_custom_target(
        cpp-format
        COMMAND ${CLANG_FORMAT_EXECUTABLE} -style=file -i ${CLANG_FORMAT_SOURCE_FILES_VALUE}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Format the C++ code with clang-format"
    )
    set_target_properties(cpp-format PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()
