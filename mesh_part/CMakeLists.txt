cmake_minimum_required(VERSION 3.16)

project(fesom_ini C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake ${CMAKE_MODULE_PATH})

# Check for exported compiler settings from main FESOM build
set(FESOM_SETTINGS_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../build/fesom_compiler_settings.cmake")
if(EXISTS ${FESOM_SETTINGS_FILE})
    message(STATUS "Using compiler settings from main FESOM build: ${FESOM_SETTINGS_FILE}")
    include(${FESOM_SETTINGS_FILE})
else()
    # Include fesom_setup when building standalone
    include(../cmake/fesom_setup.cmake)
    
    # Default compiler options when built standalone
    set(FESOM_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID})
    set(FESOM_PLATFORM_STRATEGY "notset")
    set(FESOM_ENABLE_OPENMP OFF)
    set(FESOM_ENABLE_OPENACC OFF)
endif()

# Always find MPI regardless of how we're being built
find_package(MPI REQUIRED)

# Print compiler settings being used
message(STATUS "Building mesh partitioner with compiler: ${FESOM_COMPILER_ID}")
message(STATUS "Platform strategy: ${FESOM_PLATFORM_STRATEGY}")

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}/../src)
set(sources_Fortran ${src_home}/MOD_MESH.F90 ${src_home}/oce_modules.F90 ${src_home}/gen_modules_config.F90 ${src_home}/gen_modules_partitioning.F90 ${src_home}/gen_modules_rotate_grid.F90 ${src_home}/fvom_init.F90 ${src_home}/oce_local.F90 ${src_home}/gen_comm.F90 ${src_home}/MOD_READ_BINARY_ARRAYS.F90 ${src_home}/MOD_WRITE_BINARY_ARRAYS.F90 ${src_home}/MOD_PARTIT.F90)
set(sources_C ${src_home}/fort_part.c)

# Add metis library
set(METIS_DIR ${CMAKE_CURRENT_LIST_DIR}/../lib/metis-5.1.0)
# Force METIS to build as a static library
set(SHARED FALSE CACHE BOOL "build a shared library" FORCE)
add_subdirectory(${METIS_DIR} ${PROJECT_BINARY_DIR}/metis)
include_directories(${METIS_DIR}/include)

# Build C component as a shared library
add_library(${PROJECT_NAME}_C SHARED ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(${PROJECT_NAME}_C metis)
set_target_properties(${PROJECT_NAME}_C PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Build main executable
add_executable(${PROJECT_NAME} ${sources_Fortran})

# Apply compiler flags - use exported flags if available, otherwise use preset flags
if(DEFINED FESOM_COMPILER_FLAGS)
    # Use the exact same flags as the main FESOM build
    message(STATUS "Using compiler flags from main FESOM build: ${FESOM_COMPILER_FLAGS}")
    separate_arguments(FESOM_FLAGS_LIST UNIX_COMMAND "${FESOM_COMPILER_FLAGS}")
    target_compile_options(${PROJECT_NAME} PRIVATE ${FESOM_FLAGS_LIST})
else()
    # Fall back to preset flags based on compiler ID
    message(STATUS "Using preset compiler flags for ${FESOM_COMPILER_ID}")
    
    if(${FESOM_COMPILER_ID} STREQUAL Intel OR ${FESOM_COMPILER_ID} STREQUAL IntelLLVM)
        target_compile_options(${PROJECT_NAME} PRIVATE -O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero -no-wrap-margin -fpe0)
    elseif(${FESOM_COMPILER_ID} STREQUAL GNU)
        if(${FESOM_PLATFORM_STRATEGY} STREQUAL ubuntu)
            target_compile_options(${PROJECT_NAME} PRIVATE -O2 -g -fbacktrace -ffloat-store -finit-local-zero -finline-functions -fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fallow-argument-mismatch)
        else()
            target_compile_options(${PROJECT_NAME} PRIVATE -O3 -ffloat-store -finit-local-zero -finline-functions -fimplicit-none -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none)
        endif()
        
        # Add argument mismatch flag for newer GNU compilers
        if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
            target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch)
        endif()
    elseif(${FESOM_COMPILER_ID} STREQUAL Cray)
        target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64 -N 1023 -g -G2 -O2 -hnoacc -M878)
    elseif(${FESOM_COMPILER_ID} STREQUAL NVHPC)
        target_compile_options(${PROJECT_NAME} PRIVATE -O3 -r8 -Mfree -Mallocatable=03 -Mstack_arrays -Minfo=accel -Mpreprocess)
    endif()
endif()

# Add OpenMP if enabled
if(FESOM_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_Fortran)
endif()

# Link against MPI and C library
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C MPI::MPI_Fortran)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# Make sure the executable name has .x extension for consistency
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME}.x)

# Set RPATH settings to find shared libraries
set_target_properties(${PROJECT_NAME} PROPERTIES
    INSTALL_RPATH "$ORIGIN/../lib"
    BUILD_WITH_INSTALL_RPATH TRUE
    SKIP_BUILD_RPATH FALSE
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# Installation target - install both the executable and the shared library
install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_C
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)
