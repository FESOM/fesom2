#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fesom" for configuration "Debug"
set_property(TARGET fesom APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(fesom PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_DEBUG "parms"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libfesom.so"
  IMPORTED_SONAME_DEBUG "libfesom.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS fesom )
list(APPEND _IMPORT_CHECK_FILES_FOR_fesom "${_IMPORT_PREFIX}/lib/libfesom.so" )

# Import target "parms" for configuration "Debug"
set_property(TARGET parms APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(parms PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libparms.so"
  IMPORTED_SONAME_DEBUG "libparms.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS parms )
list(APPEND _IMPORT_CHECK_FILES_FOR_parms "${_IMPORT_PREFIX}/lib/libparms.so" )

# Import target "fesom.x" for configuration "Debug"
set_property(TARGET fesom.x APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(fesom.x PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/fesom.x"
  )

list(APPEND _IMPORT_CHECK_TARGETS fesom.x )
list(APPEND _IMPORT_CHECK_FILES_FOR_fesom.x "${_IMPORT_PREFIX}/bin/fesom.x" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
