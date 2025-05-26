cmake_minimum_required(VERSION 2.8)
project(METIS)

# Set path to GKlib directory - using absolute path to avoid build directory issues
set(GKLIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")
set(SHARED FALSE CACHE BOOL "build a shared library")

#if(MSVC)
set(METIS_INSTALL FALSE)
#else()
#  set(METIS_INSTALL TRUE)
#endif()

# Configure libmetis library.
if(SHARED)
  set(METIS_LIBRARY_TYPE SHARED)
else()
  set(METIS_LIBRARY_TYPE STATIC)
endif(SHARED)

include(${GKLIB_PATH}/GKlibSystem.cmake)
# Add include directories.
include_directories(${GKLIB_PATH})
include_directories(include)
# Recursively look for CMakeLists.txt in subdirs.
add_subdirectory("include")
add_subdirectory("libmetis")
add_subdirectory("programs")
