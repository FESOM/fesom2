#===============================================================================
# tests/openmp/CMakeLists.txt - OpenMP source-level guards
#===============================================================================
#
# These tests do NOT need fesom.x to be built - they only need a configured
# tree and a Fortran compiler that can preprocess. That is deliberate: the
# defect they guard (FESOM/fesom2#444, #528) is a *compile* failure of the
# -DOPENMP_REPRODUCIBLE=ON build, so the guard has to be usable before and
# independently of building that configuration:
#
#     cmake -S . -B build -DBUILD_TESTING=ON
#     ctest --test-dir build -L openmp
#
#===============================================================================

find_package(Python3 COMPONENTS Interpreter QUIET)

if(NOT Python3_Interpreter_FOUND)
    message(WARNING
        "Python3 not found - the OpenMP ORDERED source guards "
        "(issues #444/#528) will NOT be registered.")
    return()
endif()

# Flags that make the Fortran compiler preprocess and stop.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(FORTRAN_PREPROCESS_FLAGS "-fpp -E -P")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(FORTRAN_PREPROCESS_FLAGS "-Mpreprocess -E")
else()
    set(FORTRAN_PREPROCESS_FLAGS "-cpp -E")
endif()

set(OMP_ORDERED_CHECK "${CMAKE_CURRENT_SOURCE_DIR}/check_omp_ordered.py")

# Fesom source root, derived from this file's location (tests/openmp/../..).
# Do not use CMAKE_SOURCE_DIR: in a bundled build (ecbundle) that is the
# bundle root, not the fesom tree, and src/ would not be found.
get_filename_component(FESOM_OMP_TEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)

# name -> extra arguments
function(add_omp_ordered_test TEST_NAME)
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE} "${OMP_ORDERED_CHECK}"
                --source-dir "${FESOM_OMP_TEST_ROOT}"
                --compiler "${CMAKE_Fortran_COMPILER}"
                --preprocess-flags "${FORTRAN_PREPROCESS_FLAGS}"
                ${ARGN}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        TIMEOUT 300
        LABELS "openmp;static"
    )
endfunction()

# Every !$OMP ORDERED region on the __openmp_reproducible path must sit in a
# worksharing loop that carries the ORDERED clause. This is exactly the
# compile error reported in issues #444 and #528.
add_omp_ordered_test(openmp_ordered_reproducible --mode reproducible)

# The same, for code that only the OpenIFS / IFS-interface builds compile
# (e.g. the tr_array_id==4 block in ice_fct.F90).
add_omp_ordered_test(openmp_ordered_reproducible_oifs
    --mode reproducible --extra-defines "-D__oifs -D__ifsinterface")

# The mirror invariant: neither ORDERED regions nor ORDERED clauses may leak
# into the default build, where they would serialize loops that are supposed
# to run with the OpenMP locks instead.
add_omp_ordered_test(openmp_ordered_default --mode default)

message(STATUS "Added OpenMP ORDERED source guards (issues #444/#528)")
