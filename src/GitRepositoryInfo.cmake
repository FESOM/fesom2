find_package(Git QUIET)
if(Git_FOUND)
    execute_process(
        COMMAND
            ${GIT_EXECUTABLE} rev-parse --short HEAD
        RESULT_VARIABLE
            SHORT_HASH_RESULT
        OUTPUT_VARIABLE
            SHORT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    set(SHORT_HASH "unknown")
    message("git not found, setting FESOM_GIT_SHA to: ${SHORT_HASH}")
endif()

# If running in script mode (this runs on every build)
if (CMAKE_SCRIPT_MODE_FILE)
    if (EXISTS "${SHORT_HASH_FILE}")
        file(READ ${SHORT_HASH_FILE} READ_IN_SHORT_HASH)
    else()
        set(READ_IN_SHORT_HASH "")
    endif()

    if (NOT ("${READ_IN_SHORT_HASH}" STREQUAL "${SHORT_HASH}"))
        message(STATUS "fesom git SHA is out of date")
        # This will update short_hash.txt, causing cmake to reconfigure
        file(WRITE ${SHORT_HASH_FILE} ${SHORT_HASH})
    endif()

# Else running as part of cmake configure
else()
    set(SHORT_HASH_FILE ${CMAKE_BINARY_DIR}/short_hash.txt)
    file(WRITE ${SHORT_HASH_FILE} ${SHORT_HASH})

    # The trick here is to make sure short_hash.txt is listed as a byproduct
    add_custom_target(
        git_short_hash
        BYPRODUCTS
            ${SHORT_HASH_FILE}
        COMMAND
            ${CMAKE_COMMAND}
            "-DSHORT_HASH_FILE=${SHORT_HASH_FILE}"
            "-P" "${CMAKE_CURRENT_LIST_FILE}"
        COMMENT
	"Re-determining fesom git SHA..."
        VERBATIM
        USES_TERMINAL)

    # This configure_file makes cmake reconfigure dependent on short_hash.txt
    set(FESOM_GIT_SHA ${SHORT_HASH})
    configure_file(${FESOM_ORIGINAL_VERSION_FILE} ${FESOM_GENERATED_VERSION_FILE} @ONLY)

    message(STATUS "fesom git SHA: ${SHORT_HASH}")
endif()



