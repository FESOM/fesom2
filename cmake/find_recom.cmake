# =============================================================================
# Finds the RECOM include and library directories
# 
# Inputs:
# - TOPLEVEL_DIR: defined in the top-level FESOM CMakeLists.txt
# 
# Returns:
# - RECOM_LIB_DIR: directory where librecom.a and *.mod files are found
# 
# Deniz Ural, AWI / ESM-Tools, July, 2021
# =============================================================================

# ls settings. Eg. /usr/bin/ls
# -F adds / to the directories whih is needed for regex
find_program(LS ls)
set(ls_opts "-F")

# it is assumed that RECOM is at the same level as the FESOM directory
set(search_dir "${TOPLEVEL_DIR}/../")
execute_process(COMMAND ${LS} ${ls_opts} ${search_dir}  OUTPUT_VARIABLE ls_out)
# convert multi-line output to a list for the foreach loop
string(REPLACE "\n" ";" ls_out ${ls_out})

# find RECOM library using regex. Directories with the substring 'recom' will
# be added to the candidates list
set(matches "")
foreach(item ${ls_out})
    # the directories will have a trailing slash with ls -F
    # CMake regex does not support case-insensitive search to fake it with []
    #string(REGEX MATCH ".*recom.*/$" match ${item})
    string(REGEX MATCH ".*[rR][eE][cC][oO][mM].*/$" match ${item})
    if (NOT ${match} STREQUAL "")
        list(APPEND matches ${match})
    endif()
endforeach()

# add the regex matches to the possible RECOM directories
# Since RECOM is also built with CMake, the out-of-source build will be in the 
# `build` directory
set(possible_recom_dirs "")
foreach(item ${matches})
    list(APPEND possible_recom_dirs  "${search_dir}/${item}/build/")
endforeach()

find_path(RECOM_LIB_DIR librecom.a ${possible_recom_dirs})

# check if RECOM_LIB_DIR is found. If not stop further execution.
if (${RECOM_LIB_DIR} STREQUAL RECOM_LIB_DIR-NOTFOUND)
    message(WARNING "!!! ERROR: RECOM is not found")
    message(STATUS "These matches are found for Recom")
    foreach(item ${matches})
        message(STATUS "    ${item}")
    endforeach()
    message(FATAL_ERROR "find_recom.cmake: Exiting the compilation")
else()
    message(STATUS "RECOM is found in: ${RECOM_LIB_DIR}")
endif()

