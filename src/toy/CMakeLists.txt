cmake_minimum_required(VERSION 3.16)

project(fesom_toy LANGUAGES Fortran VERSION 1.0.0)

message(STATUS "====================================================")
message(STATUS "Configuring FESOM2 Toy Ocean Model")
message(STATUS "====================================================")

# Get the toy source directory
set(toy_src_home ${CMAKE_CURRENT_SOURCE_DIR})

# Collect toy model source files (excluding the main program)
set(toy_sources_Fortran
    ${toy_src_home}/fesom_toy_module.F90
)

# Create a shared library for the toy model module
# This allows the toy model to be used independently or as a library
add_library(fesom_toy_lib SHARED ${toy_sources_Fortran})

# The toy library depends on the main FESOM library
target_link_libraries(fesom_toy_lib PUBLIC fesom)
target_link_libraries(fesom_toy_lib PRIVATE MPI::MPI_Fortran)

# Set library properties
set_target_properties(fesom_toy_lib PROPERTIES
    OUTPUT_NAME "fesom_toy"
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    LINKER_LANGUAGE Fortran
)

# Make sure module files are accessible
target_include_directories(fesom_toy_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
    $<INSTALL_INTERFACE:module/fesom_toy>
)

# Inherit compile options from parent FESOM project if available
if(TARGET fesom)
    get_target_property(FESOM_COMPILE_OPTIONS fesom COMPILE_OPTIONS)
    if(FESOM_COMPILE_OPTIONS)
        target_compile_options(fesom_toy_lib PRIVATE ${FESOM_COMPILE_OPTIONS})
    endif()

    get_target_property(FESOM_COMPILE_DEFINITIONS fesom COMPILE_DEFINITIONS)
    if(FESOM_COMPILE_DEFINITIONS)
        target_compile_definitions(fesom_toy_lib PRIVATE ${FESOM_COMPILE_DEFINITIONS})
    endif()
endif()

# Build the fesom_toy executable
add_executable(fesom_toy ${toy_src_home}/fesom_toy_main.F90)

# Link the executable against the toy library (which in turn links to fesom)
target_link_libraries(fesom_toy PUBLIC fesom_toy_lib)
target_link_libraries(fesom_toy PRIVATE MPI::MPI_Fortran)

# Installation rules
install(TARGETS fesom_toy fesom_toy_lib
        EXPORT fesom_toy-targets
        RUNTIME DESTINATION bin
            COMPONENT toy
        LIBRARY DESTINATION lib
            COMPONENT toy
        ARCHIVE DESTINATION lib
            COMPONENT toy
)

# Install Fortran module files
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
        DESTINATION module/fesom_toy
        COMPONENT toy
        FILES_MATCHING PATTERN "fesom_toy*.mod"
)

# Install the README documentation
install(FILES ${toy_src_home}/README.md
        DESTINATION share/doc/fesom_toy
        COMPONENT toy
)

# Export the toy targets for build tree
export(TARGETS fesom_toy fesom_toy_lib
       NAMESPACE fesom_toy::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-targets.cmake"
)

# Export the toy targets for install tree (for use in other CMake projects)
install(EXPORT fesom_toy-targets
        FILE fesom_toy-targets.cmake
        NAMESPACE fesom_toy::
        DESTINATION lib/cmake/fesom_toy
        COMPONENT toy
)

# Create a config file for the toy model
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/fesom_toy-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config.cmake"
    INSTALL_DESTINATION lib/cmake/fesom_toy
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# Install config files if the .in file exists
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fesom_toy-config.cmake.in")
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config-version.cmake"
        DESTINATION lib/cmake/fesom_toy
        COMPONENT toy
    )
endif()

# Display build information
message(STATUS "Toy model source files: ${toy_sources_Fortran}")
message(STATUS "Toy model library: fesom_toy_lib (shared)")
message(STATUS "Toy model executable: fesom_toy")
message(STATUS "====================================================")
