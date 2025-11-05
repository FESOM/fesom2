cmake_minimum_required(VERSION 3.16)

project(fesom_toy LANGUAGES Fortran VERSION 1.0.0)

message(STATUS "====================================================")
message(STATUS "Configuring FESOM2 Toy Ocean Model (Independent Build)")
message(STATUS "====================================================")

# Get paths
set(toy_src_home ${CMAKE_CURRENT_SOURCE_DIR})
set(fesom_src_home ${CMAKE_CURRENT_SOURCE_DIR}/../../../src)

message(STATUS "toy_src_home: ${toy_src_home}")
message(STATUS "fesom_src_home: ${fesom_src_home}")
message(STATUS "cmake path will be: ${fesom_src_home}/../cmake")

# ============================================================================
# FIND DEPENDENCIES (for standalone build)
# ============================================================================
# Check if we're building standalone or as part of main FESOM
if(NOT TARGET fesom)
    message(STATUS "Standalone build detected - finding dependencies")

    # Add cmake modules from main FESOM repository
    # Use fesom_src_home which is already set correctly for both standalone and integrated builds
    set(CMAKE_MODULE_PATH ${fesom_src_home}/../cmake ${CMAKE_MODULE_PATH})

    # Get platform strategy from environment (set by env.sh)
    if(DEFINED ENV{FESOM_PLATFORM_STRATEGY})
        set(FESOM_PLATFORM_STRATEGY $ENV{FESOM_PLATFORM_STRATEGY})
    else()
        set(FESOM_PLATFORM_STRATEGY "notset")
    endif()

    # Find MPI
    find_package(MPI REQUIRED)

    # Find NetCDF (note: FindNETCDF.cmake uses all caps)
    find_package(NETCDF REQUIRED)
    if(NetCDF_FOUND)
        set(NETCDF_C_LIBRARIES ${NetCDF_C_LIBRARIES})
        set(NETCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARIES})
        set(NETCDF_Fortran_INCLUDE_DIRECTORIES ${NetCDF_Fortran_INCLUDE_DIRS})
    endif()

    # OpenMP (optional)
    if(ENABLE_OPENMP)
        find_package(OpenMP)
    endif()
endif()

# ============================================================================
# COLLECT MINIMAL FESOM SOURCE FILES
# ============================================================================
# Only include the 78 source files actually needed by the toy model
# (excludes ice physics, coupling, biogeochemistry, etc.)

# Core Modules (MOD_*)
list(APPEND toy_fesom_sources
    ${fesom_src_home}/MOD_DYN.F90
    ${fesom_src_home}/MOD_ICE.F90
    ${fesom_src_home}/MOD_MESH.F90
    ${fesom_src_home}/MOD_PARTIT.F90
    ${fesom_src_home}/MOD_READ_BINARY_ARRAYS.F90
    ${fesom_src_home}/MOD_TRACER.F90
    ${fesom_src_home}/MOD_WRITE_BINARY_ARRAYS.F90
)

# Ocean Modules (oce_*)
list(APPEND toy_fesom_sources
    ${fesom_src_home}/oce_adv_tra_driver.F90
    ${fesom_src_home}/oce_adv_tra_fct.F90
    ${fesom_src_home}/oce_adv_tra_hor.F90
    ${fesom_src_home}/oce_adv_tra_ver.F90
    ${fesom_src_home}/oce_ale.F90
    ${fesom_src_home}/oce_ale_mixing_kpp.F90
    ${fesom_src_home}/oce_ale_mixing_pp.F90
    ${fesom_src_home}/oce_ale_pressure_bv.F90
    ${fesom_src_home}/oce_ale_ssh_splitexpl_subcycl.F90
    ${fesom_src_home}/oce_ale_tracer.F90
    ${fesom_src_home}/oce_ale_vel_rhs.F90
    ${fesom_src_home}/oce_dyn.F90
    ${fesom_src_home}/oce_fer_gm.F90
    ${fesom_src_home}/oce_mesh.F90
    ${fesom_src_home}/oce_mo_conv.F90
    ${fesom_src_home}/oce_modules.F90
    ${fesom_src_home}/oce_muscl_adv.F90
    ${fesom_src_home}/oce_setup_step.F90
    ${fesom_src_home}/oce_shortwave_pene_dbgyre.F90
    ${fesom_src_home}/oce_spp.F90
    ${fesom_src_home}/oce_tracer_mod.F90
)

# Global Modules (g_*)
list(APPEND toy_fesom_sources
    ${fesom_src_home}/cavity_param.F90
    ${fesom_src_home}/gen_events.F90
    ${fesom_src_home}/gen_forcing_init.F90
    ${fesom_src_home}/gen_halo_exchange.F90
    ${fesom_src_home}/gen_ic3d.F90
    ${fesom_src_home}/gen_interpolation.F90
    ${fesom_src_home}/gen_model_setup.F90
    ${fesom_src_home}/gen_modules_backscatter.F90
    ${fesom_src_home}/gen_modules_clock.F90
    ${fesom_src_home}/gen_modules_config.F90
    ${fesom_src_home}/gen_modules_diag.F90
    ${fesom_src_home}/gen_modules_forcing.F90
    ${fesom_src_home}/gen_modules_partitioning.F90
    ${fesom_src_home}/gen_modules_read_NetCDF.F90
    ${fesom_src_home}/gen_modules_rotate_grid.F90
    ${fesom_src_home}/gen_support.F90
)

# I/O Modules (io_*)
list(APPEND toy_fesom_sources
    ${fesom_src_home}/io_blowup.F90
    ${fesom_src_home}/io_data_strategy.F90
    ${fesom_src_home}/io_fesom_file.F90
    ${fesom_src_home}/io_gather.F90
    ${fesom_src_home}/io_meandata.F90
    ${fesom_src_home}/io_mesh_info.F90
    ${fesom_src_home}/io_netcdf_attribute_module.F90
    ${fesom_src_home}/io_netcdf_file_module.F90
    ${fesom_src_home}/io_netcdf_workaround_module.F90
    ${fesom_src_home}/io_restart.F90
    ${fesom_src_home}/io_restart_derivedtype.F90
    ${fesom_src_home}/io_restart_file_group.F90
    ${fesom_src_home}/io_scatter.F90
)

# Forcing Modules
list(APPEND toy_fesom_sources
    ${fesom_src_home}/forcing_lookahead_reader_module.F90
    ${fesom_src_home}/forcing_provider_async_module.F90
    ${fesom_src_home}/forcing_provider_netcdf_module.F90
    ${fesom_src_home}/gen_surface_forcing.F90
)

# Toy-Specific Modules
list(APPEND toy_fesom_sources
    ${fesom_src_home}/toy_channel_dbgyre.F90
    ${fesom_src_home}/toy_channel_soufflet.F90
)

# Iceberg Modules (conditionally imported, needed for compilation even though not actively used)
list(APPEND toy_fesom_sources
    ${fesom_src_home}/icb_allocate.F90
    ${fesom_src_home}/icb_coupling.F90
    ${fesom_src_home}/icb_dyn.F90
    ${fesom_src_home}/icb_elem.F90
    ${fesom_src_home}/icb_modules.F90
    ${fesom_src_home}/icb_step.F90
    ${fesom_src_home}/icb_thermo.F90
)

# Support Modules
list(APPEND toy_fesom_sources
    ${fesom_src_home}/async_threads_module.F90
    ${fesom_src_home}/command_line_options.F90
    ${fesom_src_home}/fortran_utils.F90
    ${fesom_src_home}/info_module.F90
    ${fesom_src_home}/mod_transit.F90
    ${fesom_src_home}/mpi_topology_module.F90
    ${fesom_src_home}/nvfortran_subarray_workaround.F90
    ${fesom_src_home}/solver.F90
    ${fesom_src_home}/write_step_info.F90
)

# Add version info (generated file from parent build)
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/../fesom_version_info-generated.F90")
    list(APPEND toy_fesom_sources ${CMAKE_CURRENT_BINARY_DIR}/../fesom_version_info-generated.F90)
else()
    list(APPEND toy_fesom_sources ${fesom_src_home}/fesom_version_info.F90)
endif()

# Profiler (optional)
if(FESOM_PROFILING)
    list(APPEND toy_fesom_sources ${fesom_src_home}/fesom_profiler.F90)
endif()

# CVMix modules (conditional)
# Note: We only include FESOM CVMix driver files, not CVMix library sources
# The CVMix library will be linked as an external dependency
if(CVMIX)
    # FESOM CVMix drivers
    list(APPEND toy_fesom_sources
        ${fesom_src_home}/cvmix_driver/cvmix_idemix.F90
        ${fesom_src_home}/cvmix_driver/cvmix_kinds_and_types_addon.F90
        ${fesom_src_home}/cvmix_driver/cvmix_tke.F90
        ${fesom_src_home}/cvmix_driver/cvmix_utils_addon.F90
        ${fesom_src_home}/cvmix_driver/gen_modules_cvmix_idemix.F90
        ${fesom_src_home}/cvmix_driver/gen_modules_cvmix_kpp.F90
        ${fesom_src_home}/cvmix_driver/gen_modules_cvmix_pp.F90
        ${fesom_src_home}/cvmix_driver/gen_modules_cvmix_tidal.F90
        ${fesom_src_home}/cvmix_driver/gen_modules_cvmix_tke.F90
    )
endif()

# IFS interface I/O modules (for MULTIO support, not coupling)
if(ENABLE_MULTIO)
    list(APPEND toy_fesom_sources
        ${fesom_src_home}/ifs_interface/iom.F90
        ${fesom_src_home}/ifs_interface/mpp_io.F90
    )
endif()

# Toy model module itself
list(APPEND toy_fesom_sources
    ${toy_src_home}/fesom_toy_module.F90
)

list(LENGTH toy_fesom_sources TOY_SOURCE_COUNT)
message(STATUS "Toy model will be built from ${TOY_SOURCE_COUNT} FESOM source files")

# ============================================================================
# BUILD TOY LIBRARY FROM MINIMAL SOURCES
# ============================================================================

# Create a separate module directory for the toy library to avoid conflicts
set(TOY_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/module/fesom_toy")
file(MAKE_DIRECTORY ${TOY_MODULE_DIR})

# Set Fortran module output directory for toy library compilation
set(CMAKE_Fortran_MODULE_DIRECTORY_BACKUP ${CMAKE_Fortran_MODULE_DIRECTORY})
set(CMAKE_Fortran_MODULE_DIRECTORY ${TOY_MODULE_DIR})

# Create shared library from minimal FESOM sources + toy module
add_library(fesom_toy_lib SHARED ${toy_fesom_sources})

# Restore the original module directory for other targets
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY_BACKUP})

# Add async_threads_cpp object files (from parent build)
if(TARGET async_threads_cpp)
    target_sources(fesom_toy_lib PRIVATE $<TARGET_OBJECTS:async_threads_cpp>)
else()
    # If building standalone, compile async_threads_cpp with -fPIC for shared library
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    add_subdirectory(${fesom_src_home}/async_threads_cpp ${CMAKE_CURRENT_BINARY_DIR}/async_threads_cpp)
    target_sources(fesom_toy_lib PRIVATE $<TARGET_OBJECTS:async_threads_cpp>)
endif()

# Set library properties
set_target_properties(fesom_toy_lib PROPERTIES
    OUTPUT_NAME "fesom_toy"
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    LINKER_LANGUAGE Fortran
    Fortran_MODULE_DIRECTORY ${TOY_MODULE_DIR}
)

# Include directories - use toy's own module directory
target_include_directories(fesom_toy_lib PUBLIC
    $<BUILD_INTERFACE:${TOY_MODULE_DIR}>
    $<INSTALL_INTERFACE:module/fesom_toy>
)

# Also need to search the toy module directory for dependencies
target_include_directories(fesom_toy_lib PRIVATE ${TOY_MODULE_DIR})

# Link against external dependencies (inherited from parent or found locally)
target_link_libraries(fesom_toy_lib PRIVATE MPI::MPI_Fortran)
target_link_libraries(fesom_toy_lib PRIVATE ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES})

# Link against CVMix if enabled
if(CVMIX AND TARGET cvmix_external)
    target_link_libraries(fesom_toy_lib PRIVATE cvmix_external)
endif()

# Link against MULTIO if enabled
if(ENABLE_MULTIO AND TARGET multio-fapi)
    target_link_libraries(fesom_toy_lib PRIVATE multio-fapi)
    target_compile_definitions(fesom_toy_lib PRIVATE __MULTIO)
endif()

# OpenMP support
if(ENABLE_OPENMP AND TARGET OpenMP::OpenMP_Fortran)
    target_link_libraries(fesom_toy_lib PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Set compiler flags (inherit from fesom or set directly for standalone)
if(TARGET fesom)
    # Building as part of main FESOM - inherit flags
    get_target_property(FESOM_COMPILE_OPTIONS fesom COMPILE_OPTIONS)
    if(FESOM_COMPILE_OPTIONS)
        target_compile_options(fesom_toy_lib PRIVATE ${FESOM_COMPILE_OPTIONS})
    endif()

    get_target_property(FESOM_COMPILE_DEFINITIONS fesom COMPILE_DEFINITIONS)
    if(FESOM_COMPILE_DEFINITIONS)
        target_compile_definitions(fesom_toy_lib PRIVATE ${FESOM_COMPILE_DEFINITIONS})
    endif()
else()
    # Standalone build - use same compiler flags as main FESOM (from src/CMakeLists.txt)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        if(FESOM_PLATFORM_STRATEGY STREQUAL "ubuntu")
            message(STATUS "Allowing type mismatches on Ubuntu for CI Testing")
            target_compile_options(fesom_toy_lib PRIVATE -O2 -g -fbacktrace -ffloat-store -finit-local-zero -finline-functions -fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fallow-argument-mismatch -cpp)
        else()
            target_compile_options(fesom_toy_lib PRIVATE -O3 -ffloat-store -finit-local-zero -finline-functions -fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp)
        endif()
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
            target_compile_options(fesom_toy_lib PRIVATE -fallow-argument-mismatch)
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        target_compile_options(fesom_toy_lib PRIVATE -O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero -no-wrap-margin -fpe0 -fpp)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
            target_compile_options(fesom_toy_lib PRIVATE -no-prec-sqrt -ip)
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
        target_compile_options(fesom_toy_lib PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O2 -hnoacc -M878)
    endif()
endif()

# Additional compile definitions
if(CVMIX)
    target_compile_definitions(fesom_toy_lib PRIVATE __cvmix)
endif()
if(FESOM_PROFILING)
    target_compile_definitions(fesom_toy_lib PRIVATE FESOM_PROFILING)
endif()
if(DISABLE_MULTITHREADING)
    target_compile_definitions(fesom_toy_lib PRIVATE DISABLE_MULTITHREADING)
endif()
if(VERBOSE)
    target_compile_definitions(fesom_toy_lib PRIVATE VERBOSE)
endif()
if(OPENMP_REPRODUCIBLE)
    target_compile_definitions(fesom_toy_lib PRIVATE __openmp_reproducible)
endif()

# Include NetCDF directories
target_include_directories(fesom_toy_lib PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES})

# ============================================================================
# BUILD TOY EXECUTABLE
# ============================================================================

add_executable(fesom_toy ${toy_src_home}/fesom_toy_main.F90)
target_link_libraries(fesom_toy PUBLIC fesom_toy_lib)
target_link_libraries(fesom_toy PRIVATE MPI::MPI_Fortran)
# Link against C++ standard library (needed for async_threads_cpp)
target_link_libraries(fesom_toy PRIVATE stdc++)

# ============================================================================
# INSTALLATION
# ============================================================================

install(TARGETS fesom_toy fesom_toy_lib
        EXPORT fesom_toy-targets
        RUNTIME DESTINATION bin
            COMPONENT toy
        LIBRARY DESTINATION lib
            COMPONENT toy
        ARCHIVE DESTINATION lib
            COMPONENT toy
)

# Install Fortran module files (from toy's own module directory)
install(DIRECTORY ${TOY_MODULE_DIR}/
        DESTINATION module/fesom_toy
        COMPONENT toy
        FILES_MATCHING PATTERN "*.mod"
)

# Install documentation
install(FILES ${toy_src_home}/../README.md
        DESTINATION share/doc/fesom_toy
        COMPONENT toy
)

# Export for build tree
export(TARGETS fesom_toy fesom_toy_lib
       NAMESPACE fesom_toy::
       FILE "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-targets.cmake"
)

# Export for install tree
install(EXPORT fesom_toy-targets
        FILE fesom_toy-targets.cmake
        NAMESPACE fesom_toy::
        DESTINATION lib/cmake/fesom_toy
        COMPONENT toy
)

# Package configuration
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

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fesom_toy-config.cmake.in")
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/fesom_toy-config-version.cmake"
        DESTINATION lib/cmake/fesom_toy
        COMPONENT toy
    )
endif()

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "====================================================")
message(STATUS "Toy Model Configuration Summary:")
message(STATUS "  Source files: ${TOY_SOURCE_COUNT}")
message(STATUS "  Library: libfesom_toy.so (independent, minimal FESOM)")
message(STATUS "  Executable: fesom_toy")
message(STATUS "  CVMIX: ${CVMIX}")
message(STATUS "  MULTIO: ${ENABLE_MULTIO}")
message(STATUS "  OpenMP: ${ENABLE_OPENMP}")
message(STATUS "  Profiling: ${FESOM_PROFILING}")
message(STATUS "====================================================")
