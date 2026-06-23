#===============================================================================
# CheckFesomRun.cmake - robust pass/fail verification for FESOM CTest runs
#===============================================================================
#
# Provides check_fesom_run(), used by the generated per-test run scripts
# (see cmake/FesomTesting.cmake) and by the mesh partitioner test.
#
# WHY THIS EXISTS
# ---------------
# The previous harness decided pass/fail from the process exit code alone, and
# even accepted exit code 1 as success. That is unsound for FESOM because:
#   * par_ex(..., abort) aborts via MPI_ABORT(COMM, 1) -> exit code 1.
#   * FESOM frequently prints "ERROR ..." and then a bare Fortran `stop`, which
#     exits 0 under gfortran; par_ex without `abort` also returns 0.
# So a crashed or aborted run could be reported green.
#
# A run is considered PASSED only if ALL of these hold:
#   (1) exit code == 0
#   (2) every SUCCESS_MARKER (if any) is present in the captured output
#   (3) no FAILURE_SIGNATURE is present in the captured output
#   (4) every REQUIRED_ARTIFACT path exists
# Otherwise the function raises FATAL_ERROR (-> CTest failure) and dumps a
# diagnostic summary plus the tail of the captured output.
#
# Usage (from a `cmake -P` script):
#   include("/abs/path/to/cmake/CheckFesomRun.cmake")
#   check_fesom_run(
#       NAME               my_test
#       RESULT             "${test_result}"
#       OUTPUT_LOG         "${run_dir}/test_output.log"
#       ERROR_LOG          "${run_dir}/test_error.log"
#       SUCCESS_MARKERS    "fesom should stop with exit status = 0"
#       REQUIRED_ARTIFACTS "${result_dir}/sst.fesom.1948.nc"
#   )
#===============================================================================

if(COMMAND check_fesom_run)
    return()
endif()

function(check_fesom_run)
    set(oneValueArgs NAME RESULT OUTPUT_LOG ERROR_LOG)
    set(multiValueArgs SUCCESS_MARKERS FAILURE_SIGNATURES REQUIRED_ARTIFACTS)
    cmake_parse_arguments(CR "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT CR_NAME)
        set(CR_NAME "fesom_run")
    endif()

    # Default failure signatures shared by FESOM standalone runs and the tools.
    # Caller may override by passing FAILURE_SIGNATURES (pass an empty string to
    # disable signature scanning entirely).
    if(NOT DEFINED CR_FAILURE_SIGNATURES)
        set(CR_FAILURE_SIGNATURES
            "Run finished unexpectedly!"   # par_ex abort path (rank 0)
            "ERROR:"                       # FESOM convention: write(*,*) 'ERROR: ...'
            "MPI_ABORT"
            "MPI_Abort"
            "forrtl:"                      # Intel runtime errors
            "Backtrace for this error"     # gfortran crash backtrace
            "Segmentation fault"
            "segmentation fault")
    endif()

    # Gather combined stdout+stderr from the captured log files.
    set(_combined "")
    foreach(_log "${CR_OUTPUT_LOG}" "${CR_ERROR_LOG}")
        if(_log AND EXISTS "${_log}")
            file(READ "${_log}" _content)
            string(APPEND _combined "${_content}\n")
        endif()
    endforeach()

    set(_failures "")

    # (1) exit code must be exactly 0
    if(NOT "${CR_RESULT}" STREQUAL "0")
        list(APPEND _failures "non-zero exit code: '${CR_RESULT}'")
    endif()

    # (2) every success marker must be present
    foreach(_marker IN LISTS CR_SUCCESS_MARKERS)
        if(NOT "${_marker}" STREQUAL "")
            string(FIND "${_combined}" "${_marker}" _pos)
            if(_pos EQUAL -1)
                list(APPEND _failures "missing success marker: '${_marker}'")
            endif()
        endif()
    endforeach()

    # (3) no failure signature may appear
    foreach(_sig IN LISTS CR_FAILURE_SIGNATURES)
        if(NOT "${_sig}" STREQUAL "")
            string(FIND "${_combined}" "${_sig}" _pos)
            if(NOT _pos EQUAL -1)
                list(APPEND _failures "failure signature found in output: '${_sig}'")
            endif()
        endif()
    endforeach()

    # (4) every required artifact must exist
    foreach(_art IN LISTS CR_REQUIRED_ARTIFACTS)
        if(NOT "${_art}" STREQUAL "")
            if(NOT EXISTS "${_art}")
                list(APPEND _failures "missing required artifact: ${_art}")
            endif()
        endif()
    endforeach()

    if(_failures)
        list(LENGTH _failures _nfail)
        message("")
        message("================ FESOM test FAILED: ${CR_NAME} ================")
        foreach(_f IN LISTS _failures)
            message("  - ${_f}")
        endforeach()
        if(CR_OUTPUT_LOG)
            message("  stdout log: ${CR_OUTPUT_LOG}")
        endif()
        if(CR_ERROR_LOG)
            message("  stderr log: ${CR_ERROR_LOG}")
        endif()
        # Tail of captured output to aid diagnosis without opening the log files.
        string(LENGTH "${_combined}" _len)
        if(_len GREATER 2000)
            math(EXPR _start "${_len} - 2000")
            string(SUBSTRING "${_combined}" ${_start} -1 _tail)
            message("  --- last 2000 chars of output ---\n${_tail}")
        elseif(_len GREATER 0)
            message("  --- captured output ---\n${_combined}")
        endif()
        message(FATAL_ERROR "Test ${CR_NAME} FAILED (${_nfail} check(s) failed)")
    else()
        message(STATUS "Test ${CR_NAME} PASSED (exit 0; success markers present; no failure signatures; required artifacts present)")
    endif()
endfunction()
