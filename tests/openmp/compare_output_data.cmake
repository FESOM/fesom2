#===============================================================================
# compare_output_data.cmake - bit-for-bit comparison of two FESOM result dirs
#===============================================================================
#
# cmake -DNCDUMP=<ncdump> -DDIR_A=<results> -DDIR_B=<results> -DFILES="a.nc;b.nc" -P compare_output_data.cmake
#
# For every file, the data section of `ncdump -p 9,17` (17 significant digits,
# enough to round-trip a double exactly) must be identical in both directories.
# The header is skipped: it carries run-specific global attributes.
#===============================================================================

foreach(_var NCDUMP DIR_A DIR_B FILES)
    if(NOT DEFINED ${_var})
        message(FATAL_ERROR "compare_output_data.cmake: ${_var} not set")
    endif()
endforeach()

set(_differ "")
foreach(_f IN LISTS FILES)
    foreach(_side A B)
        set(_path "${DIR_${_side}}/${_f}")
        if(NOT EXISTS "${_path}")
            message(FATAL_ERROR "missing output file: ${_path}")
        endif()
        execute_process(COMMAND "${NCDUMP}" -p 9,17 "${_path}"
                        OUTPUT_VARIABLE _dump RESULT_VARIABLE _rc)
        if(NOT _rc EQUAL 0)
            message(FATAL_ERROR "ncdump failed on ${_path}")
        endif()
        string(FIND "${_dump}" "\ndata:" _pos)
        string(SUBSTRING "${_dump}" ${_pos} -1 _data_${_side})
    endforeach()
    if(_data_A STREQUAL _data_B)
        message(STATUS "identical: ${_f}")
    else()
        list(APPEND _differ "${_f}")
    endif()
endforeach()

if(_differ)
    message(FATAL_ERROR "outputs differ between\n  ${DIR_A}\n  ${DIR_B}\nin: ${_differ}")
endif()
