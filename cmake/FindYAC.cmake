# Check for conflict between BUILD_YAC and yac_DIR environment variable
if(BUILD_YAC AND DEFINED ENV{yac_DIR})
    message(FATAL_ERROR
        "Cannot use both BUILD_YAC=ON and yac_DIR environment variable.\n"
        "Please choose one approach:\n"
        "  - Use BUILD_YAC=ON to automatically build YAC from source, OR\n"
        "  - Set yac_DIR environment variable to use existing YAC installation")
endif()

# If BUILD_YAC is enabled, BuildYAC.cmake will set all necessary variables
if(BUILD_YAC)
    return()
endif()

IF( DEFINED ENV{yac_DIR} )
  SET( yac_DIR "$ENV{yac_DIR}" )
ENDIF()

set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${yac_DIR}/lib/pkgconfig:${yac_DIR}/src/pkgconfig")

# The Fortran driver uses the module "yac" (src/mci/yac_module.F90 in the
# YAC sources) together with yac_fread_config_yaml, yac_fdef_calendar,
# yac_fget_comp_comm and yac_fenddef. All of them require YAC 3.2.0 or
# newer, which is also the release that added yac-mci.pc and deprecated
# yac.pc.
set(YAC_MIN_VERSION 3.2.0)

find_package(PkgConfig QUIET)

# newer versions of yac dont expose "yac" but "yac-mci"
pkg_check_modules(PC_yac QUIET "yac-mci >= ${YAC_MIN_VERSION}")
if(NOT PC_yac_FOUND)
  pkg_check_modules(PC_yac QUIET "yac >= ${YAC_MIN_VERSION}")
endif()

if(NOT PC_yac_FOUND)
  # diagnose the failure: no pkg-config, no YAC, or a YAC too old
  if(NOT PKG_CONFIG_FOUND)
    message(FATAL_ERROR
        "FESOM_COUPLING=yac needs pkg-config to locate YAC "
        "(or use -DBUILD_YAC=ON to build YAC from source)")
  endif()
  pkg_check_modules(PC_yac_unversioned QUIET yac-mci)
  if(NOT PC_yac_unversioned_FOUND)
    pkg_check_modules(PC_yac_unversioned QUIET yac)
  endif()
  if(PC_yac_unversioned_FOUND)
    message(FATAL_ERROR
        "FESOM_COUPLING=yac requires YAC >= ${YAC_MIN_VERSION}, "
        "but found YAC ${PC_yac_unversioned_VERSION} in ${yac_DIR}")
  else()
    message(FATAL_ERROR
        "FESOM_COUPLING=yac requires YAC >= ${YAC_MIN_VERSION}, but no "
        "yac-mci.pc or yac.pc was found in PKG_CONFIG_PATH.\n"
        "Please choose one approach:\n"
        "  - Set the yac_DIR environment variable to a YAC installation, OR\n"
        "  - Use BUILD_YAC=ON to build YAC from source")
  endif()
endif()

message(STATUS "Found YAC ${PC_yac_VERSION}")

find_path(YAC_Fortran_INCLUDE_DIRECTORIES NAMES yac.mod
          HINTS ${PC_yac_INCLUDE_DIRS} ${yac_DIR}/src/mci)

# pkg-config hands us YAC as raw -L/-l flags, from which
# CMAKE_INSTALL_RPATH_USE_LINK_PATH cannot derive an RPATH entry. Export
# the library directories so the caller can add them explicitly, which a
# YAC installation with shared libraries needs both to link and to run.
set(YAC_LIBRARY_DIRECTORIES ${PC_yac_LIBRARY_DIRS})

find_library(YAC_LIBRARY yac yac_mci HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/src/mci ${yac_DIR}/src)
#find_library(YAC_CLAPACK_LIBRARY yac_clapack HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/clapack)
#find_library(YAC_MTIME_LIBRARY yac_mtime HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/mtime)

#list(REMOVE_ITEM PC_yac_LINK_LIBRARIES yac yac_clapack yac_mtime)

#message(FATAL_ERROR "${PC_yac_LDFLAGS}")
set(YAC_Fortran_LIBRARIES "-L${yac_DIR}/src" "-L${yac_DIR}/src/mci" ${PC_yac_LDFLAGS})
