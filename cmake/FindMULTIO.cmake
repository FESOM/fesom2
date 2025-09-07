# FindMULTIO.cmake

include(FindPackageHandleStandardArgs)

# Use the environment variable as a hint
set(MULTIO_HINT_PATH $ENV{MULTIO_INSTALL_PATH})

# Try to find the library
find_library(MULTIO_FAPI_LIBRARY
    NAMES multio-fapi # Adjust this if the library name is different
    HINTS ${MULTIO_HINT_PATH}/lib
)

# Try to find the dependency library
find_library(MULTIO_API_LIBRARY
    NAMES multio-api
    HINTS ${MULTIO_HINT_PATH}/lib
)

# Try to find the dependency library
find_library(MULTIO_LIBRARY
    NAMES multio
    HINTS ${MULTIO_HINT_PATH}/lib
)

# Try to find the Fortran module path
find_path(MULTIO_MODULE_PATH
    NAMES multio_api.mod # Replace <module_name> with an actual module name you expect to find
    HINTS ${MULTIO_HINT_PATH}/module ${MULTIO_HINT_PATH}/multio/module
)


# Aggregate the libraries for easier linking
set(MULTIO_LIBRARIES ${MULTIO_FAPI_LIBRARY} ${MULTIO_API_LIBRARY} ${MULTIO_LIBRARY})


# Handle the results
find_package_handle_standard_args(MULTIO 
    REQUIRED_VARS MULTIO_LIBRARIES MULTIO_MODULE_PATH
    FOUND_VAR MULTIO_FOUND
)

# If found, set the MULTIO_LIBRARIES and MULTIO_INCLUDE_DIRS variables for easy use
if(MULTIO_FOUND)
    set(MULTIO_INCLUDE_DIRS ${MULTIO_MODULE_PATH})
endif()

# Mark variables as advanced
mark_as_advanced(MULTIO_LIBRARY MULTIO_MODULE_PATH)

