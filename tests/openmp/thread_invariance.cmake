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

# Compare the runs with A and B threads; an optional third argument names a solver
# variant (run directories openmp_pi_mpi2_<variant>_omp<N>).
function(_omp_add_compare A B)
    set(_v "")
    if(ARGC GREATER 2)
        set(_v "${ARGV2}_")
    endif()
    set(_name openmp_pi_mpi2_${_v}identical_omp${A}_omp${B})
    add_test(NAME ${_name}
        COMMAND ${CMAKE_COMMAND}
            -DNCDUMP=${NCDUMP_EXECUTABLE}
            -DDIR_A=${CMAKE_CURRENT_BINARY_DIR}/openmp_pi_mpi2_${_v}omp${A}/results
            -DDIR_B=${CMAKE_CURRENT_BINARY_DIR}/openmp_pi_mpi2_${_v}omp${B}/results
            "-DFILES=${_omp_files}"
            -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_output_data.cmake)
    set_tests_properties(${_name} PROPERTIES
        LABELS "openmp"
        TIMEOUT 300
        FIXTURES_REQUIRED "openmp_pi_mpi2_${_v}omp${A};openmp_pi_mpi2_${_v}omp${B}")
endfunction()

_omp_add_compare(1 4)
_omp_add_compare(1 8)
_omp_add_compare(8 8_rerun)

# Set KEY = VALUE in one namelist of a configured run directory.
function(_omp_set_namelist RUN_DIR FILE KEY VALUE)
    file(READ "${RUN_DIR}/${FILE}" _nl)
    string(REGEX REPLACE "(\n[ \t]*${KEY}[ \t]*=)[^!\n]*" "\\1 ${VALUE} " _nl "${_nl}")
    file(WRITE "${RUN_DIR}/${FILE}" "${_nl}")
endfunction()

# Solver variants that the default configuration does not reach, each run with 1
# and 4 threads per rank and compared bit for bit:
#   mevp  modified EVP sea-ice rheology (whichEVP=1)
#   aevp  adaptive EVP sea-ice rheology (whichEVP=2)
#   se    split-explicit barotropic subcycling (use_ssh_se_subcycl)
set(_omp_variants mevp aevp se)
set(_omp_variant_mevp namelist.ice whichEVP 1)
set(_omp_variant_aevp namelist.ice whichEVP 2)
set(_omp_variant_se   namelist.dyn use_ssh_se_subcycl .true.)

foreach(_v IN LISTS _omp_variants)
    foreach(_threads 1 4)
        set(_name openmp_pi_mpi2_${_v}_omp${_threads})
        add_fesom_test_with_options(${_name}
            "pi" "96" "${FESOM_OPENMP_TEST_DAYS}" "d" "${FESOM_OPENMP_TEST_DAYS}" "d" "96" ".true." ".false."
            MPI_TEST
            NP 2
            OMP_THREADS ${_threads}
            LABEL openmp
            TIMEOUT ${_omp_timeout}
        )
        _omp_set_output("${CMAKE_CURRENT_BINARY_DIR}/${_name}")
        _omp_set_namelist("${CMAKE_CURRENT_BINARY_DIR}/${_name}" ${_omp_variant_${_v}})
        set_tests_properties(${_name} PROPERTIES FIXTURES_SETUP ${_name})
    endforeach()
    _omp_add_compare(1 4 ${_v})
endforeach()

message(STATUS "Added OpenMP thread-invariance tests (pi, 2 ranks x 1/4/8 threads, ${FESOM_OPENMP_TEST_DAYS} d; variants: ${_omp_variants})")
