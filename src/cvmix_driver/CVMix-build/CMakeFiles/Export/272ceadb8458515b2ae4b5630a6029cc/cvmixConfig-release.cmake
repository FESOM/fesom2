#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cvmix_static" for configuration "Release"
set_property(TARGET cvmix_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cvmix_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcvmix.a"
  )

list(APPEND _cmake_import_check_targets cvmix_static )
list(APPEND _cmake_import_check_files_for_cvmix_static "${_IMPORT_PREFIX}/lib/libcvmix.a" )

# Import target "cvmix_shared" for configuration "Release"
set_property(TARGET cvmix_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cvmix_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcvmix.so"
  IMPORTED_SONAME_RELEASE "libcvmix.so"
  )

list(APPEND _cmake_import_check_targets cvmix_shared )
list(APPEND _cmake_import_check_files_for_cvmix_shared "${_IMPORT_PREFIX}/lib/libcvmix.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
