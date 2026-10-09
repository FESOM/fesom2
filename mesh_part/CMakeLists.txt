# ----------------------------------------------------------------------
# fesom_ini — FESOM2's mesh partitioner (mesh_part subproject), ported and made
# mixed-mesh capable (triangle + quad elements, use_depthfile='nemo', the
# appended dist_N sections; see docs/HANDOFF.md "Hybrid (mixed triangle + quad)
# meshes").
#
# Reads namelist.config (paths, geometry, run_config, machine) from the working
# directory and writes edges.out, edge_tri.out, edgenum.out, elvls.out,
# elvls_raw.out, nlvls.out and dist_N/ (incl. dist_N/format.out) INTO the mesh
# directory (MeshPath).
#
# The sources under src/ derive from FESOM2's (/home/a/a270029/fesom2.0/src). They define
# their own MOD_MESH / MOD_PARTIT / o_PARAM ..., so they get a private module
# directory and never see (or clash with) FESOM3's modules. Compile flags and
# the C definitions of fort_part.c are FESOM2's mesh_part/CMakeLists.txt.
# ----------------------------------------------------------------------
project(fesom_ini C Fortran)

set(src_home   ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(METIS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../lib/metis-5.1.0)

# Private module directory (FESOM2's mod_mesh/mod_partit vs FESOM3's).
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

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

if(NOT IS_DIRECTORY ${METIS_ROOT})
    message(FATAL_ERROR "METIS not found at ${METIS_ROOT}.")
endif()
message(STATUS "fesom_ini: METIS from ${METIS_ROOT}")
# EXCLUDE_FROM_ALL: only libmetis is built (as a dependency), not METIS's programs.
add_subdirectory(${METIS_ROOT} ${CMAKE_CURRENT_BINARY_DIR}/metis EXCLUDE_FROM_ALL)

# C library for mesh partitioning (FESOM2 definitions, incl. the METIS seed).
add_library(fesom_ini_C STATIC ${sources_C})
target_include_directories(fesom_ini_C PRIVATE ${METIS_ROOT}/include)
target_compile_definitions(fesom_ini_C PRIVATE
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
    METIS_VERSION=5
    PART_WEIGHTED
    METISRANDOMSEED=35243
)
target_link_libraries(fesom_ini_C metis)

add_executable(fesom_ini ${sources_Fortran})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(fesom_ini PRIVATE
        -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip
    )
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fesom_ini PRIVATE
        -fdefault-real-8 -ffree-line-length-none
    )
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        target_compile_options(fesom_ini PRIVATE -fallow-argument-mismatch)
    endif()
endif()
target_link_libraries(fesom_ini PUBLIC fesom_ini_C MPI::MPI_Fortran)
set_target_properties(fesom_ini PROPERTIES LINKER_LANGUAGE Fortran)
