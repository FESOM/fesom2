#===============================================================================
# thread_invariance.cmake - OpenMP thread-count invariance on the pi mesh
#===============================================================================
#
# The same pi run (2 MPI ranks) with 1, 4 and 8 OpenMP threads per rank, and the
# 8-thread run repeated. All outputs must be bit-identical: the model's threaded
# loops accumulate in a mesh-defined order, so neither the thread count nor the
# thread scheduling may change the answer.
#
#     ctest -L openmp
#
# FESOM_OPENMP_TEST_DAYS sets the run length. The default (1 day) is the pull-request
# check: a scheduling-dependent sum already changes the first day at every node. The
# manual workflow fesom2_openmp_invariance.yml runs longer (e.g. 10 or 30 days).
#
# Needs an OpenMP build (ENABLE_OPENMP=ON), MPI tests and ncdump.
#===============================================================================

if(NOT ENABLE_OPENMP OR NOT ENABLE_MPI_TESTS OR FESOM_COUPLED OR OIFS_COUPLED OR USE_YAC)
    return()
endif()

find_program(NCDUMP_EXECUTABLE ncdump)
if(NOT NCDUMP_EXECUTABLE)
    message(WARNING "ncdump not found - the OpenMP thread-invariance tests will NOT be registered.")
    return()
endif()

set(FESOM_OPENMP_TEST_DAYS 1 CACHE STRING "Length in days of the OpenMP thread-invariance runs")
math(EXPR _omp_timeout "300 + 120 * ${FESOM_OPENMP_TEST_DAYS}")

set(_omp_fields sst sss ssh uice vice a_ice m_ice temp salt u v w)
set(_omp_files "")
foreach(_f IN LISTS _omp_fields)
    list(APPEND _omp_files "${_f}.fesom.1948.nc")
endforeach()

# Daily means of the listed fields in double precision, so that a difference in
# the last bit is visible.
function(_omp_set_output RUN_DIR)
    set(_list "")
    foreach(_f IN LISTS _omp_fields)
        if(_list)
            string(APPEND _list ",\n           ")
        endif()
        string(APPEND _list "'${_f}',1, 'd', 8")
    endforeach()
    file(READ "${RUN_DIR}/namelist.io" _io)
    string(REGEX REPLACE "io_list[ \t]*=[^/]*/" "io_list = ${_list}\n/" _io "${_io}")
    file(WRITE "${RUN_DIR}/namelist.io" "${_io}")
endfunction()

set(_omp_runs 1 4 8 8_rerun)
foreach(_run IN LISTS _omp_runs)
    string(REGEX REPLACE "_rerun$" "" _threads "${_run}")
    set(_name openmp_pi_mpi2_omp${_run})
    add_fesom_test_with_options(${_name}
        "pi" "96" "${FESOM_OPENMP_TEST_DAYS}" "d" "${FESOM_OPENMP_TEST_DAYS}" "d" "96" ".true." ".false."
        MPI_TEST
        NP 2
        OMP_THREADS ${_threads}
        LABEL openmp
        TIMEOUT ${_omp_timeout}
    )
    _omp_set_output("${CMAKE_CURRENT_BINARY_DIR}/${_name}")
    set_tests_properties(${_name} PROPERTIES FIXTURES_SETUP ${_name})
endforeach()

function(_omp_add_compare A B)
    set(_name openmp_pi_mpi2_identical_omp${A}_omp${B})
    add_test(NAME ${_name}
        COMMAND ${CMAKE_COMMAND}
            -DNCDUMP=${NCDUMP_EXECUTABLE}
            -DDIR_A=${CMAKE_CURRENT_BINARY_DIR}/openmp_pi_mpi2_omp${A}/results
            -DDIR_B=${CMAKE_CURRENT_BINARY_DIR}/openmp_pi_mpi2_omp${B}/results
            "-DFILES=${_omp_files}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_output_data.cmake)
    set_tests_properties(${_name} PROPERTIES
        LABELS "openmp"
        TIMEOUT 300
        FIXTURES_REQUIRED "openmp_pi_mpi2_omp${A};openmp_pi_mpi2_omp${B}")
endfunction()

_omp_add_compare(1 4)
_omp_add_compare(1 8)
_omp_add_compare(8 8_rerun)

message(STATUS "Added OpenMP thread-invariance tests (pi, 2 ranks x 1/4/8 threads, ${FESOM_OPENMP_TEST_DAYS} d)")
