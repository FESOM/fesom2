# To a large degree this follows the ecbuild setup
message("[${PROJECT_NAME}]")

### Setup project

if(CMAKE_VERSION VERSION_LESS 3.21)
  get_property(not_top DIRECTORY PROPERTY PARENT_DIRECTORY)
  if(NOT not_top)
    set(PROJECT_IS_TOP_LEVEL true)
  endif()
endif()

set(BUILD_SHARED_LIBS ON CACHE BOOL "Default to using shared libs")
set(CMAKE_LINK_DEPENDS_NO_SHARED ON) # relink of downstream libraries not required when shared library is rebuilt

# Default build type: Release
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")

# Set location to look for find_package modules
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/../cmake ${CMAKE_MODULE_PATH})

# Set Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/${PROJECT_NAME} )

# Build-dir destinations
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set install directories
include(GNUInstallDirs)
if( NOT INSTALL_BIN_DIR )
   set( INSTALL_BIN_DIR ${CMAKE_INSTALL_BINDIR} )
endif()
if( NOT INSTALL_LIB_DIR )
   set( INSTALL_LIB_DIR ${CMAKE_INSTALL_LIBDIR} )
endif()
if( NOT INSTALL_INCLUDE_DIR )
   set( INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR} )
endif()
set( INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/${PROJECT_NAME} )

mark_as_advanced( INSTALL_BIN_DIR )
mark_as_advanced( INSTALL_LIB_DIR )
mark_as_advanced( INSTALL_INCLUDE_DIR )
mark_as_advanced( INSTALL_CMAKE_DIR )

# make sure nothing breaks if INSTALL_LIB_DIR is not lib
if( NOT INSTALL_LIB_DIR STREQUAL "lib" AND NOT EXISTS ${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR} )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink lib ${CMAKE_BINARY_DIR}/${INSTALL_LIB_DIR} )
endif()

# for macosx use @rpath in a target’s install name (CMP0042)
set( CMAKE_MACOSX_RPATH ON )

# add the automatic parts to RPATH which point to dirs outside build tree
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH   TRUE  )

# use RPATHs for the build tree
set( CMAKE_SKIP_BUILD_RPATH              FALSE )

# build with *relative* rpaths by default
option( ENABLE_RPATHS           "when installing insert RPATHS into binaries"     ON  )
option( ENABLE_RELATIVE_RPATHS  "try to use relative RPATHS, including build dir" ON  )
if( ENABLE_RELATIVE_RPATHS )
    set( CMAKE_BUILD_WITH_INSTALL_RPATH  TRUE )
else()
    # in case the RPATH is absolute, the install RPATH cannot be set
    # at build-time since it breaks the build tree dynamic links
    set( CMAKE_BUILD_WITH_INSTALL_RPATH  FALSE )
endif()


foreach( p LIB BIN INCLUDE DATA CMAKE )
   set( var INSTALL_${p}_DIR )
   set( ${PROJECT_NAME}_FULL_INSTALL_${p}_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}_DIR}"
       CACHE INTERNAL "${PROJECT_NAME} ${p} full install path" )
endforeach()


function( _path_append var path )
  list( FIND ${var} ${path} _found )
  if( _found EQUAL "-1" )
    list( APPEND ${var} ${path})
  endif()
  set( ${var} "${${var}}" PARENT_SCOPE ) #
endfunction()

function( _make_relative_rpath_entry entry var )

   if( CMAKE_SYSTEM_NAME MATCHES "Darwin" )
      set( ${var} "@loader_path/${entry}" PARENT_SCOPE )

   elseif( CMAKE_SYSTEM_NAME MATCHES "FreeBSD|Linux|SunOS" )
      set( ${var} "$ORIGIN/${entry}" PARENT_SCOPE )

   elseif( CMAKE_SYSTEM_NAME MATCHES "AIX" ) # always relative to executable path
      set( ${var} "${entry}" PARENT_SCOPE )

   else()
      set( ${var} "${CMAKE_INSTALL_PREFIX}/${entry}" PARENT_SCOPE )

   endif()
endfunction()

macro( append_to_rpath RPATH_DIRS )

   foreach( RPATH_DIR ${RPATH_DIRS} )

        if( NOT ${RPATH_DIR} STREQUAL "" )

            file( TO_CMAKE_PATH ${RPATH_DIR} RPATH_DIR ) # sanitize the path

            if( IS_ABSOLUTE ${RPATH_DIR} )
               _path_append( CMAKE_INSTALL_RPATH "${RPATH_DIR}" )
            else()
                _make_relative_rpath_entry( "${RPATH_DIR}" rpath_dir_rel )
                _path_append( CMAKE_INSTALL_RPATH ${rpath_dir_rel} )
            endif()

     endif()

   endforeach()

endmacro( append_to_rpath )


if( ENABLE_RPATHS )
   if( ENABLE_RELATIVE_RPATHS )
      file( RELATIVE_PATH relative_rpath ${CMAKE_INSTALL_PREFIX}/${INSTALL_BIN_DIR} ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR} )
      append_to_rpath( ${relative_rpath} )
   else() # make rpaths absolute
      append_to_rpath( "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}" )
   endif()
endif()

# put the include dirs which are in the source or build tree
# before all other include dirs, so the headers in the sources
# are prefered over the already installed ones (since cmake 2.4.1)
set( CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON )

# Set -fPIC flag etc.
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

include(fesom_export)
