find_package(Git QUIET)
if(Git_FOUND)
   execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                   OUTPUT_VARIABLE FESOM_GIT_SHA
                   OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
   set(FESOM_GIT_SHA "unknown")
   message("git not found, setting FESOM_GIT_SHA to: ${FESOM_GIT_SHA}")
endif()

configure_file(${FESOM_ORIGINAL_VERSION_FILE} ${FESOM_GENERATED_VERSION_FILE} @ONLY)
