# AddPreprocessTarget.cmake
# CMake module to add a target for preprocessing Fortran source files

# Function to add a preprocess target to a CMake project
# Usage: add_preprocess_target(TARGET_NAME SOURCE_FILES [EXCLUDE_PATTERN pattern])
function(add_preprocess_target TARGET_NAME SOURCE_FILES)
    # Parse additional arguments
    set(options "")
    set(oneValueArgs EXCLUDE_PATTERN)
    set(multiValueArgs "")
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Create a list of source files to preprocess
    set(preprocess_sources ${SOURCE_FILES})
    
    # If an exclude pattern is provided, filter out matching files
    if(DEFINED ARG_EXCLUDE_PATTERN)
        list(FILTER preprocess_sources EXCLUDE REGEX "${ARG_EXCLUDE_PATTERN}")
    endif()

    # Add custom target for preprocessing
    add_custom_target(${TARGET_NAME}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/preprocessed
        COMMAND ${CMAKE_COMMAND} -E echo "Preprocessing Fortran sources to ${CMAKE_BINARY_DIR}/preprocessed"
        COMMENT "Preprocessing Fortran source files"
    )

    # Get all compile definitions from the main target
    get_target_property(TARGET_COMPILE_DEFS ${PROJECT_NAME} COMPILE_DEFINITIONS)

    # Get all compile options from the main target
    get_target_property(TARGET_COMPILE_OPTIONS ${PROJECT_NAME} COMPILE_OPTIONS)

    # Add custom command for each Fortran source file
    foreach(source_file ${preprocess_sources})
        get_filename_component(fname ${source_file} NAME)
        
        # Build the command with proper arguments based on compiler type
        if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
            set(preprocess_cmd ${CMAKE_Fortran_COMPILER} -cpp -E)
        elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM")
            set(preprocess_cmd ${CMAKE_Fortran_COMPILER} -fpp -E)
        elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
            set(preprocess_cmd ${CMAKE_Fortran_COMPILER} -Mpreprocess -E)
        else()
            # Default to GNU syntax
            set(preprocess_cmd ${CMAKE_Fortran_COMPILER} -cpp -E)
        endif()
        
        # Add compile definitions
        if(TARGET_COMPILE_DEFS)
            foreach(def ${TARGET_COMPILE_DEFS})
                list(APPEND preprocess_cmd -D${def})
            endforeach()
        endif()
        
        # Add compile options (filtering out ones that might cause issues)
        if(TARGET_COMPILE_OPTIONS)
            foreach(opt ${TARGET_COMPILE_OPTIONS})
                if(NOT opt MATCHES "-O[0-3]" AND NOT opt MATCHES "-g" AND NOT opt MATCHES "-march" AND NOT opt MATCHES "-mtune")
                    list(APPEND preprocess_cmd ${opt})
                endif()
            endforeach()
        endif()
        
        # Add source file and output file
        list(APPEND preprocess_cmd ${source_file} -o ${CMAKE_BINARY_DIR}/preprocessed/${fname})
        
        add_custom_command(
            TARGET ${TARGET_NAME}
            COMMAND ${preprocess_cmd}
            DEPENDS ${source_file}
        )
    endforeach()

endfunction()
