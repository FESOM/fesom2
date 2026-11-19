#===============================================================================
# tests/mpi_reproducibility/CMakeLists.txt - MPI rank-count reproducibility (pi)
#===============================================================================
#
# The pi mesh run on 2 and 8 ranks (dist_2, dist_8) must give
#   * a bit-identical initial state (temp_M1/salt_M1 in the restart after one
#     step hold the state the first step started from),
#   * restarts and daily means after FESOM_MPI_REPRO_DAYS within
#     FESOM_MPI_REPRO_RTOL of each other (relative to the field maximum). The
#     EVP sea-ice stresses sigma11/12/22 are left out: the subcycled rheology
#     amplifies round-off to a few percent within a day.
# With BUILD_MESHPARTITIONER the same checks run on partitions freshly made by
# fesom_meshpart, and the fresh 8-rank initial state is compared with the
# committed one.
#
# Label mpi_repro holds the checks above. Label mpi_bitwise holds the strict
# bit-for-bit comparison of the full runs, which fails until all global
# reductions are independent of the decomposition.
#
#     cmake -DBUILD_TESTING=ON -DENABLE_MPI_REPRO_TESTS=ON [-DBUILD_MESHPARTITIONER=ON]
#     ctest -L mpi_repro
#     ctest -L mpi_bitwise
#
# Needs MPI tests and a Python 3 with netCDF4.
#===============================================================================

if(NOT ENABLE_MPI_REPRO_TESTS)
    return()
endif()
if(NOT ENABLE_MPI_TESTS OR FESOM_COUPLED OR OIFS_COUPLED OR USE_YAC)
    message(STATUS "MPI reproducibility tests need a standalone build with MPI tests - not registered")
    return()
endif()

find_package(Python3 COMPONENTS Interpreter QUIET)
set(_repro_python "")
if(Python3_Interpreter_FOUND)
    foreach(_flags "" "-s")
        execute_process(COMMAND ${Python3_EXECUTABLE} ${_flags} -c "import netCDF4, numpy"
                        RESULT_VARIABLE _nc_ok OUTPUT_QUIET ERROR_QUIET)
        if(_nc_ok EQUAL 0)
            set(_repro_python ${Python3_EXECUTABLE} ${_flags})
            break()
        endif()
    endforeach()
endif()
if(NOT _repro_python)
    message(WARNING "No Python 3 with netCDF4 - the MPI reproducibility tests will NOT be registered.")
    return()
endif()

set(FESOM_MPI_REPRO_DAYS 1 CACHE STRING "Length in days of the MPI reproducibility runs")
set(FESOM_MPI_REPRO_RTOL 1e-3 CACHE STRING "Relative tolerance of the 2 vs 8 rank comparison after FESOM_MPI_REPRO_DAYS")
math(EXPR _repro_timeout "300 + 120 * ${FESOM_MPI_REPRO_DAYS}")

set(_repro_fields sst sss ssh uice vice a_ice m_ice temp salt u v w)
set(_repro_ic_vars "temp_M1,salt_M1")
set(_repro_close_exclude "sigma11,sigma12,sigma22")

function(_repro_set_output RUN_DIR UNIT)
    set(_list "")
    foreach(_f IN LISTS _repro_fields)
        if(_list)
            string(APPEND _list ",\n           ")
        endif()
        string(APPEND _list "'${_f}',1, '${UNIT}', 8")
    endforeach()
    file(READ "${RUN_DIR}/namelist.io" _io)
    string(REGEX REPLACE "io_list[ \t]*=[^/]*/" "io_list = ${_list}\n/" _io "${_io}")
    file(WRITE "${RUN_DIR}/namelist.io" "${_io}")
endfunction()

# _repro_run(<name> <mesh> <np> step|day [DATA_DIR <dir>] [FIXTURE <partition fixture>])
function(_repro_run NAME MESH NP KIND)
    cmake_parse_arguments(R "" "DATA_DIR;FIXTURE" "" ${ARGN})
    set(_extra "")
    if(R_DATA_DIR)
        set(_extra DATA_DIR "${R_DATA_DIR}")
    endif()
    if(KIND STREQUAL "step")
        set(_len 1)
        set(_unit s)
    else()
        set(_len ${FESOM_MPI_REPRO_DAYS})
        set(_unit d)
    endif()
    add_fesom_test_with_options(${NAME}
        "${MESH}" "96" "${_len}" "${_unit}" "${_len}" "${_unit}" "96" ".true." ".false."
        MPI_TEST
        FRESH_START
        NP ${NP}
        LABEL mpi_repro
        TIMEOUT ${_repro_timeout}
        ${_extra}
    )
    _repro_set_output("${CMAKE_CURRENT_BINARY_DIR}/${NAME}" ${_unit})
    set_tests_properties(${NAME} PROPERTIES FIXTURES_SETUP ${NAME})
    if(R_FIXTURE)
        set_tests_properties(${NAME} PROPERTIES FIXTURES_REQUIRED ${R_FIXTURE})
    endif()
endfunction()

# _repro_compare(<name> <label> <ref run> <test run> <args to compare_nc_close.py>...)
function(_repro_compare NAME LABEL REF TEST)
    add_test(NAME ${NAME}
        COMMAND ${_repro_python} "${FESOM_TESTING_ROOT}/tools/compare_nc_close.py"
                "${CMAKE_CURRENT_BINARY_DIR}/${REF}/results"
                "${CMAKE_CURRENT_BINARY_DIR}/${TEST}/results"
                ${ARGN})
    set_tests_properties(${NAME} PROPERTIES
        LABELS "${LABEL}"
        TIMEOUT 300
        FIXTURES_REQUIRED "${REF};${TEST}")
endfunction()

function(_repro_compare_exact NAME REF TEST)
    add_test(NAME ${NAME}
        COMMAND ${_repro_python} "${FESOM_TESTING_ROOT}/tools/compare_nc_exact.py"
                "${CMAKE_CURRENT_BINARY_DIR}/${REF}/results"
                "${CMAKE_CURRENT_BINARY_DIR}/${TEST}/results")
    set_tests_properties(${NAME} PROPERTIES
        LABELS "mpi_bitwise"
        TIMEOUT 300
        FIXTURES_REQUIRED "${REF};${TEST}")
endfunction()

# committed partitions
foreach(_np 2 8)
    _repro_run(mpi_repro_pi_step_np${_np} pi ${_np} step)
    _repro_run(mpi_repro_pi_day_np${_np}  pi ${_np} day)
endforeach()
_repro_compare(mpi_repro_pi_ic_np2_np8 mpi_repro
    mpi_repro_pi_step_np2 mpi_repro_pi_step_np8 --vars ${_repro_ic_vars})
_repro_compare(mpi_repro_pi_day_close_np2_np8 mpi_repro
    mpi_repro_pi_day_np2 mpi_repro_pi_day_np8 --rtol ${FESOM_MPI_REPRO_RTOL} --exclude ${_repro_close_exclude})
_repro_compare_exact(mpi_bitwise_pi_day_np2_np8 mpi_repro_pi_day_np2 mpi_repro_pi_day_np8)

# partitions made by fesom_meshpart
if(BUILD_MESHPARTITIONER)
    set(_fresh "${CMAKE_BINARY_DIR}/tests/mpi_repro_data")
    file(MAKE_DIRECTORY "${_fresh}/MESHES")
    file(GLOB _entries LIST_DIRECTORIES true "${FESOM_TESTING_ROOT}/tests/data/*")
    foreach(_e IN LISTS _entries)
        get_filename_component(_n "${_e}" NAME)
        if(NOT _n STREQUAL "MESHES" AND IS_DIRECTORY "${_e}" AND NOT EXISTS "${_fresh}/${_n}")
            file(CREATE_LINK "${_e}" "${_fresh}/${_n}" SYMBOLIC)
        endif()
    endforeach()

    foreach(_np 2 8)
        set(_mesh pi_mpirepro_${_np})
        set(_workdir "${CMAKE_BINARY_DIR}/tests/mpi_reproducibility/meshes/${_mesh}")
        file(MAKE_DIRECTORY "${_workdir}")
        if(NOT EXISTS "${_fresh}/MESHES/${_mesh}")
            file(CREATE_LINK "${_workdir}" "${_fresh}/MESHES/${_mesh}" SYMBOLIC)
        endif()
        add_test(NAME mpi_repro_partition_pi_np${_np}
            COMMAND ${CMAKE_COMMAND}
                -DSOURCE_DIR=${FESOM_TESTING_ROOT}
                -DBUILD_DIR=${CMAKE_BINARY_DIR}
                -DMESH_NAME=pi
                -DNUM_PROCESSES=${_np}
                -DMESH_WORKDIR=${_workdir}
                -P ${FESOM_TESTING_ROOT}/tests/integration/mesh_partition.cmake)
        set_tests_properties(mpi_repro_partition_pi_np${_np} PROPERTIES
            LABELS mpi_repro
            TIMEOUT 300
            FIXTURES_SETUP mpi_repro_partition_pi_np${_np})

        _repro_run(mpi_repro_pi_fresh_step_np${_np} ${_mesh} ${_np} step
            DATA_DIR "${_fresh}" FIXTURE mpi_repro_partition_pi_np${_np})
        _repro_run(mpi_repro_pi_fresh_day_np${_np}  ${_mesh} ${_np} day
            DATA_DIR "${_fresh}" FIXTURE mpi_repro_partition_pi_np${_np})
    endforeach()

    _repro_compare(mpi_repro_pi_fresh_ic_np2_np8 mpi_repro
        mpi_repro_pi_fresh_step_np2 mpi_repro_pi_fresh_step_np8 --vars ${_repro_ic_vars})
    _repro_compare(mpi_repro_pi_ic_committed_fresh_np8 mpi_repro
        mpi_repro_pi_step_np8 mpi_repro_pi_fresh_step_np8 --vars ${_repro_ic_vars})
    _repro_compare(mpi_repro_pi_fresh_day_close_np2_np8 mpi_repro
        mpi_repro_pi_fresh_day_np2 mpi_repro_pi_fresh_day_np8 --rtol ${FESOM_MPI_REPRO_RTOL} --exclude ${_repro_close_exclude})
    _repro_compare_exact(mpi_bitwise_pi_fresh_day_np2_np8 mpi_repro_pi_fresh_day_np2 mpi_repro_pi_fresh_day_np8)
endif()

message(STATUS "Added MPI reproducibility tests (pi, 2 vs 8 ranks, ${FESOM_MPI_REPRO_DAYS} d, rtol ${FESOM_MPI_REPRO_RTOL})")
