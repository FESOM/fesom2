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

set(SHORT_HASH_FILE ${CMAKE_BINARY_DIR}/short_hash.txt)
if (EXISTS "${SHORT_HASH_FILE}")
    file(READ ${SHORT_HASH_FILE} READ_IN_SHORT_HASH)
else()
    set(READ_IN_SHORT_HASH "")
endif()

message(STATUS "FESOM_GIT_SHA: '${FESOM_GIT_SHA}'")
if (NOT ("${READ_IN_SHORT_HASH}" STREQUAL "${FESOM_GIT_SHA}"))
    message(STATUS "previous FESOM_GIT_SHA: '${READ_IN_SHORT_HASH}' regenerating ${FESOM_GENERATED_VERSION_FILE}")
    # This will update short_hash.txt, causing cmake to reconfigure
    file(WRITE ${SHORT_HASH_FILE} ${FESOM_GIT_SHA})
    configure_file(${FESOM_ORIGINAL_VERSION_FILE} ${FESOM_GENERATED_VERSION_FILE} @ONLY)
endif()

