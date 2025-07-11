cmake_minimum_required(VERSION 3.16)
project(fesom_meshpart C Fortran)

# Find required packages
find_package(MPI REQUIRED)

# Determine source directory - works for both standalone and integrated builds
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # This is a standalone build
    set(src_home ${CMAKE_SOURCE_DIR}/../src)
    set(METIS_ROOT ${CMAKE_SOURCE_DIR}/../lib/metis-5.1.0)
    message(STATUS "Building as standalone mesh partitioner")
    
    # Set default install prefix if not specified
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/.. CACHE PATH "Install path prefix" FORCE)
    endif()
    
    # Add uninstall target
    if(NOT TARGET uninstall)
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cmake_uninstall.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
            IMMEDIATE @ONLY)
        add_custom_target(uninstall
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    endif()
else()
    # This is part of the main FESOM build
    set(src_home ${CMAKE_SOURCE_DIR}/src)
    set(METIS_ROOT ${CMAKE_SOURCE_DIR}/lib/metis-5.1.0)
    message(STATUS "Building as part of FESOM main build")
endif()

message(STATUS "Source directory: ${src_home}")
message(STATUS "METIS root: ${METIS_ROOT}")

# Source files for the mesh tools
set(sources_Fortran 
    ${src_home}/MOD_MESH.F90 
    ${src_home}/oce_modules.F90 
    ${src_home}/gen_modules_config.F90 
    ${src_home}/gen_modules_partitioning.F90 
    ${src_home}/gen_modules_rotate_grid.F90 
    ${src_home}/fvom_init.F90 
    ${src_home}/oce_local.F90 
    ${src_home}/gen_comm.F90 
    ${src_home}/MOD_READ_BINARY_ARRAYS.F90 
    ${src_home}/MOD_WRITE_BINARY_ARRAYS.F90 
    ${src_home}/MOD_PARTIT.F90
)
set(sources_C ${src_home}/fort_part.c)

# Check if metis exists
if(EXISTS ${METIS_ROOT} AND IS_DIRECTORY ${METIS_ROOT})
    message(STATUS "Using METIS from: ${METIS_ROOT}")
    add_subdirectory(${METIS_ROOT} ${PROJECT_BINARY_DIR}/metis)
    include_directories(${METIS_ROOT}/include)
else()
    message(FATAL_ERROR "METIS not found at ${METIS_ROOT}. Please ensure the METIS library is available.")
endif()

# C library for mesh partitioning
add_library(${PROJECT_NAME}_C ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE
    USE_MPI 
    REAL=double
    DBL 
    HAS_BLAS 
    FORTRAN_UNDERSCORE 
    VOID_POINTER_SIZE_8 
    SGI 
    LINUX 
    UNDER_ 
    MPI2
)
target_link_libraries(${PROJECT_NAME}_C metis)

# Main executable
add_executable(${PROJECT_NAME} ${sources_Fortran})

# Add install target for standalone build
if(NOT DEFINED FESOM_SOURCE_DIR)
    install(TARGETS ${PROJECT_NAME} 
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            COMPONENT tools)
    
    # Install the C library as well
    install(TARGETS ${PROJECT_NAME}_C
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib
            COMPONENT tools)
endif()

# Set compiler flags based on the compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(${PROJECT_NAME} PRIVATE 
        -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip
    )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${PROJECT_NAME} PRIVATE 
        -fdefault-real-8 -ffree-line-length-none
    )
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch)
    endif()
endif()
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C MPI::MPI_Fortran)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

# Install the executable
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    COMPONENT tools
)
