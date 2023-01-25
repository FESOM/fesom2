IF( DEFINED ENV{yac_DIR} )
  SET( yac_DIR "$ENV{yac_DIR}" )
ENDIF()

set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${yac_DIR}/lib/pkgconfig:${yac_DIR}/src/pkgconfig")

find_package(PkgConfig QUIET)
PKG_CHECK_MODULES(PC_yac yac)

find_path(YAC_Fortran_INCLUDE_DIRECTORIES mo_yac_finterface.mod HINTS ${PC_yac_INCLUDE_DIRS} ${yac_DIR}/src)

find_library(YAC_LIBRARY yac HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/src)
find_library(YAC_CLAPACK_LIBRARY yac_clapack HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/clapack)
find_library(YAC_MTIME_LIBRARY yac_mtime HINTS ${PC_yac_LINK_LIBRARIES} ${yac_DIR}/mtime)

list(REMOVE_ITEM PC_yac_LINK_LIBRARIES yac yac_clapack yac_mtime)

set(YAC_Fortran_LIBRARIES ${YAC_LIBRARY} ${YAC_CLAPACK_LIBRARY} ${YAC_MTIME_LIBRARY} ${PC_yac_LINK_LIBRARIES} ${PC_yac_LDFLAGS_OTHER})
