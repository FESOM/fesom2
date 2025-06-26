cmake_minimum_required(VERSION 3.16)
project(fesom_meshpart C Fortran)

# Find required packages
find_package(MPI REQUIRED)

# Get our source files from the main src directory
set(src_home ${CMAKE_SOURCE_DIR}/src)

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

# Add metis library - using the one from the main source tree
add_subdirectory(${CMAKE_SOURCE_DIR}/lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
include_directories(${CMAKE_SOURCE_DIR}/lib/metis-5.1.0/include)

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