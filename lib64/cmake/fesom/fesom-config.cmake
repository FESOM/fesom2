# This file will get included upon `find_package(fesom ... )`


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was fesom-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

### insert definitions for IMPORTED targets
if(NOT fesom_BINARY_DIR)
  find_file(fesom_TARGETS_FILE
    NAMES fesom-targets.cmake
    HINTS ${CMAKE_CURRENT_LIST_DIR} 
    NO_DEFAULT_PATH)
  if(fesom_TARGETS_FILE)
    include(${fesom_TARGETS_FILE})
  endif()
endif()

include(CMakeFindDependencyMacro)

if( NOT (TARGET MPI::MPI_C AND TARGET MPI::MPI_Fortran) )
  enable_language(C)
  enable_language(Fortran)
  find_dependency(MPI COMPONENTS C Fortran)
endif()

if( OFF )
  if( NOT TARGET OpenMP::OpenMP_Fortran )
    enable_language(Fortran)
    find_dependency(OpenMP COMPONENTS Fortran)
  endif()
  set( fesom_HAVE_OPENMP 1 )
else()
  set( fesom_HAVE_OPENMP 0 )
endif()

if( OFF )
  set( fesom_HAVE_IFS_INTERFACE 1 )
else()
  set( fesom_HAVE_IFS_INTERFACE 0 )
endif()

if( OFF )
  find_dependency( multio   HINTS ${CMAKE_CURRENT_LIST_DIR}/../multio    )
  set( fesom_HAVE_MULTIO 1 )
else()
  set( fesom_HAVE_MULTIO 0 )
endif()
