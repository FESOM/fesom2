cmake_minimum_required(VERSION 3.16)
project(fesom3 Fortran)

# ----------------------------------------------------------------------
# FESOM3 — clean-architecture Fortran reimplementation of FESOM2.
# v1 anchor build (bit-identity oracle = FESOM2 v2.7.3): Intel + double
# precision. See cmake/fesom_flags.cmake and docs/LESSONS.md.
# ----------------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Default to the anchor (Release) build unless told otherwise.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Build type (Release = anchor; Debug = +checks, FP-identical)" FORCE)
endif()

# Precision options (mutually exclusive; DP is default = anchor).
option(USE_SINGLE_PRECISION "Build with single precision (WP=4)" OFF)
option(USE_HALF_PRECISION   "Build with half precision (WP=2, FP16; NVHPC only)" OFF)
if(USE_HALF_PRECISION AND USE_SINGLE_PRECISION)
    message(FATAL_ERROR "USE_HALF_PRECISION and USE_SINGLE_PRECISION are mutually exclusive.")
endif()
if(USE_HALF_PRECISION)
    message(STATUS "FESOM3 precision: HALF (WP=2)")
elseif(USE_SINGLE_PRECISION)
    message(STATUS "FESOM3 precision: SINGLE (WP=4)")
else()
    message(STATUS "FESOM3 precision: DOUBLE (WP=8, anchor)")
endif()

include(fesom_flags)

find_package(MPI REQUIRED COMPONENTS Fortran)

# ----------------------------------------------------------------------
# netCDF (Fortran) — needed from M2.10 (forcing read). Discover via nf-config
# (on PATH from env/levante module load). This MUST be the same spack
# netcdf-fortran the FESOM2 oracle links (4.5.3) so the raw real(4) slices read
# byte-identically. Exposed as an INTERFACE target `fesom_netcdf`.
# ----------------------------------------------------------------------
find_program(NF_CONFIG nf-config)
if(NOT NF_CONFIG)
    message(FATAL_ERROR "nf-config not found on PATH — load the netcdf-fortran module (env/levante.dkrz.de/shell.*).")
endif()
execute_process(COMMAND ${NF_CONFIG} --includedir
                OUTPUT_VARIABLE NETCDF_F_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${NF_CONFIG} --flibs
                OUTPUT_VARIABLE NETCDF_F_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(NETCDF_F_FLIBS_LIST UNIX_COMMAND "${NETCDF_F_FLIBS}")
message(STATUS "FESOM3 netCDF-Fortran include: ${NETCDF_F_INCLUDE}")
message(STATUS "FESOM3 netCDF-Fortran flibs  : ${NETCDF_F_FLIBS}")
# Bake the netcdf-fortran + netcdf-c lib dirs into the binaries' RPATH so they run
# without LD_LIBRARY_PATH (the spack modules set compile paths but not runtime).
execute_process(COMMAND ${NF_CONFIG} --prefix OUTPUT_VARIABLE NF_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
find_program(NC_CONFIG nc-config)
set(NETCDF_RPATH "${NF_PREFIX}/lib")
if(NC_CONFIG)
    execute_process(COMMAND ${NC_CONFIG} --prefix OUTPUT_VARIABLE NC_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    list(APPEND NETCDF_RPATH "${NC_PREFIX}/lib")
endif()
add_library(fesom_netcdf INTERFACE)
target_include_directories(fesom_netcdf INTERFACE ${NETCDF_F_INCLUDE})
target_link_libraries(fesom_netcdf INTERFACE ${NETCDF_F_FLIBS_LIST})
foreach(_rp ${NETCDF_RPATH})
    target_link_options(fesom_netcdf INTERFACE "-Wl,-rpath,${_rp}")
endforeach()

# ----------------------------------------------------------------------
# liblz4 (OPTIONAL) — the first Zarr output compressor (M9, Task 0.3). We call LZ4 via Fortran
# iso_c_binding (our own interfaces, no C header), so only the LINK library matters. The runtime
# soname liblz4.so.1 is a system lib (/lib64); the dev symlink liblz4.so lives in spack. Discover
# it; if absent, build none-codec-only (mod_io_zarr error-stops on an lz4 request — HAVE_LZ4 guard).
# Exposed as the INTERFACE target `fesom_lz4`, mirroring `fesom_netcdf`.
# ----------------------------------------------------------------------
find_library(LZ4_LIB NAMES lz4
    HINTS
        /sw/spack-levante/lz4-1.9.4-fhvgx6/lib
        /sw/spack-levante/lz4-1.9.4-jjx6m3/lib
        /sw/spack-levante/lz4-1.9.4-htgmak/lib
        /usr/lib64 /lib64)
if(LZ4_LIB)
    message(STATUS "FESOM3 lz4: ${LZ4_LIB} (HAVE_LZ4)")
    add_library(fesom_lz4 INTERFACE)
    target_link_libraries(fesom_lz4 INTERFACE ${LZ4_LIB})
    get_filename_component(_lz4_dir ${LZ4_LIB} DIRECTORY)
    target_link_options(fesom_lz4 INTERFACE "-Wl,-rpath,${_lz4_dir}")
else()
    message(STATUS "FESOM3 lz4: NOT FOUND — building none-codec-only (no HAVE_LZ4)")
endif()

# Put .mod files in a private dir so the scanner resolves USE across subdirs.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ----------------------------------------------------------------------
# Library: everything under src/ except the driver entry points.
# GLOB + CONFIGURE_DEPENDS so newly added module files are picked up on
# rebuild without re-running cmake by hand. CMake's Fortran dependency
# scanner orders the modules (params -> types -> infra -> ...) automatically.
# ----------------------------------------------------------------------
set(FESOM3_LIB_DIRS params types infra mesh oce ice forcing io step)
set(fesom3_lib_sources "")
foreach(_d ${FESOM3_LIB_DIRS})
    file(GLOB _s CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/${_d}/*.F90)
    list(APPEND fesom3_lib_sources ${_s})
endforeach()

add_library(fesom3 STATIC ${fesom3_lib_sources})
apply_fesom_compile_flags(fesom3)
target_link_libraries(fesom3 PUBLIC MPI::MPI_Fortran fesom_netcdf)
if(LZ4_LIB)
    target_link_libraries(fesom3 PUBLIC fesom_lz4)
    target_compile_definitions(fesom3 PRIVATE HAVE_LZ4)
endif()
set_target_properties(fesom3 PROPERTIES LINKER_LANGUAGE Fortran)

# ----------------------------------------------------------------------
# Driver executables: one per file in src/drivers/.
# ----------------------------------------------------------------------
file(GLOB fesom3_drivers CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/*.F90)
foreach(_drv ${fesom3_drivers})
    get_filename_component(_name ${_drv} NAME_WE)
    add_executable(${_name} ${_drv})
    apply_fesom_compile_flags(${_name})
    target_link_libraries(${_name} PUBLIC fesom3)
    set_target_properties(${_name} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()

# ----------------------------------------------------------------------
# fesom_ini: FESOM2's mesh partitioner (mesh_part/, METIS in lib/metis-5.1.0).
# Writes edges/levels and dist_N/ into the mesh directory; see mesh_part/CMakeLists.txt.
# ----------------------------------------------------------------------
add_subdirectory(mesh_part)

# ----------------------------------------------------------------------
# Tests: one executable per test/test_*.F90, registered with ctest.
# add_fesom_test(<name> <nranks>) builds + runs under mpirun.
# ----------------------------------------------------------------------
enable_testing()

# Oversubscribe-safe launcher (login-node runs of small rank counts).
set(FESOM3_MPIEXEC ${MPIEXEC_EXECUTABLE})
if(NOT FESOM3_MPIEXEC)
    set(FESOM3_MPIEXEC mpirun)
endif()
# The DKRZ UCX runtime settings (env/levante shells) target compute-node
# InfiniBand and fail on login nodes. Self-tests are tiny single-node runs, so
# default to the ob1/vader (shared-memory) path. Override for batch if needed.
# vader's single-copy mechanism (knem, the Levante default) is switched off: Open MPI
# 4.1.2 caches knem registrations by virtual address and misses glibc heap trims (brk),
# so a send buffer re-allocated at the same address after a trim is read from the OLD
# pinned pages -- right count, payload of an earlier message (LESSONS L61; standalone
# evidence: tools/run_mpi_knem_stale_repro.sh). env.sh sets the same for all runs.
set(FESOM3_MPI_LAUNCH_FLAGS "--mca;pml;ob1;--mca;btl;self,vader;--mca;btl_vader_single_copy_mechanism;none;--oversubscribe"
    CACHE STRING "Extra mpirun flags for ctest self-tests")
# a build tree configured before the single-copy switch keeps its cached vader flags:
# add the switch there (a cache without vader, e.g. a UCX override, is left alone)
if(FESOM3_MPI_LAUNCH_FLAGS MATCHES "vader" AND NOT FESOM3_MPI_LAUNCH_FLAGS MATCHES "single_copy_mechanism")
    list(INSERT FESOM3_MPI_LAUNCH_FLAGS 0 --mca btl_vader_single_copy_mechanism none)
    set(FESOM3_MPI_LAUNCH_FLAGS "${FESOM3_MPI_LAUNCH_FLAGS}" CACHE STRING
        "Extra mpirun flags for ctest self-tests" FORCE)
    message(STATUS "FESOM3_MPI_LAUNCH_FLAGS: added btl_vader_single_copy_mechanism none")
endif()

function(add_fesom_test test_name nranks)
    # Build the executable once; allow registering the same test at several rank
    # counts (add_fesom_test called multiple times for one test_name).
    if(NOT TARGET ${test_name})
        add_executable(${test_name} ${CMAKE_SOURCE_DIR}/test/${test_name}.F90)
        apply_fesom_compile_flags(${test_name})
        target_link_libraries(${test_name} PUBLIC fesom3)
        if(TARGET fesom3_test_util)
            target_link_libraries(${test_name} PUBLIC fesom3_test_util)
        endif()
        set_target_properties(${test_name} PROPERTIES
            LINKER_LANGUAGE Fortran
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    endif()
    add_test(NAME ${test_name}_np${nranks}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${nranks}
                     $<TARGET_FILE:${test_name}>)
endfunction()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt)
    add_subdirectory(test)
endif()
