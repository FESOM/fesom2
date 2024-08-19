# NETCDF_Fortran_INCLUDE_DIRECTORIES
# NETCDF_Fortran_LIBRARIES
# NETCDF_C_INCLUDE_DIRECTORIES
# NETCDF_C_LIBRARIES


# include file netcdf.inc and search symbol e.g. NCCRE for Fortran 77 NetCDF
# include file netcdf.mod and search symbol e.g. nf90_create for Fortran 90 NetCDF (this symbol is not found on hlrn with intel ftn)
if(CMAKE_Fortran_COMPILER_LOADED)
	include(CheckFortranFunctionExists)
	check_fortran_function_exists(NCCRE HAVE_Fortran_NETCDF)
	
	if(HAVE_Fortran_NETCDF)
		set(NETCDF_Fortran_INCLUDE_DIRECTORIES "")
		set(NETCDF_Fortran_LIBRARIES "")
		set(NETCDF_Fortran_FOUND 1)
	else()
		find_path(NETCDF_Fortran_INCLUDE_DIRECTORIES netcdf.inc HINTS $ENV{NETCDF_ROOT}/include $ENV{NETCDF_DIR}/include $ENV{NETCDF4_DIR}/include ENV NETCDF_Fortran_INCLUDE_DIRECTORIES)
		find_library(NETCDF_Fortran_LIBRARIES netcdff HINTS ${NETCDF_Fortran_INCLUDE_DIRECTORIES}/../lib)
		if( NETCDF_Fortran_INCLUDE_DIRECTORIES AND NETCDF_Fortran_LIBRARIES )
			set(NETCDF_Fortran_FOUND 1)
		endif()
	endif()	
endif()

# if one wants to link to the static NetCDF, also link these libraries for C: hdf5, hdf5_hl, curl
if(CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
	include(CheckFunctionExists)
	check_function_exists(nc__create HAVE_C_NETCDF)

	if(HAVE_C_NETCDF)
		set(NETCDF_C_INCLUDE_DIRECTORIES "")
		set(NETCDF_C_LIBRARIES "")
		set(NETCDF_C_FOUND 1)
	else()
		find_path(NETCDF_C_INCLUDE_DIRECTORIES netcdf.h HINTS $ENV{NETCDF_ROOT}/include $ENV{NETCDF_DIR}/include $ENV{NETCDF4_DIR}/include ENV NETCDF_C_INCLUDE_DIRECTORIES)
		find_library(NETCDF_C_LIBRARIES netcdf HINTS ${NETCDF_C_INCLUDE_DIRECTORIES}/../lib)
		if( NETCDF_C_INCLUDE_DIRECTORIES AND NETCDF_C_LIBRARIES )
			set(NETCDF_C_FOUND 1)
		endif()
	endif()
endif()

if(CMAKE_CXX_COMPILER_LOADED)
	find_path(NETCDF_CXX_INCLUDE_DIRECTORIES netcdf HINTS $ENV{NETCDF_ROOT}/include $ENV{NETCDF_DIR}/include $ENV{NETCDF4_DIR} ENV NETCDF_CXX_INCLUDE_DIRECTORIES)
	# the cray toolchain (e.g. hlrn) disables dynamic linking by default. to enable it at build time do e.g. "CRAYPE_LINK_TYPE=dynamic make".	
	find_library(NETCDF_CXX_LIBRARIES NAMES netcdf_c++4 netcdf-cxx4 HINTS ${NETCDF_CXX_INCLUDE_DIRECTORIES}/../lib)
	if(NETCDF_CXX_INCLUDE_DIRECTORIES AND NETCDF_C_INCLUDE_DIRECTORIES)
		list(APPEND NETCDF_CXX_INCLUDE_DIRECTORIES ${NETCDF_C_INCLUDE_DIRECTORIES})
	endif()
	if(NETCDF_CXX_LIBRARIES AND NETCDF_C_LIBRARIES)
		list(APPEND NETCDF_CXX_LIBRARIES ${NETCDF_C_LIBRARIES})
	endif()
	if( NETCDF_CXX_INCLUDE_DIRECTORIES AND NETCDF_CXX_LIBRARIES )
		set(NETCDF_CXX_FOUND 1)
	endif()
endif()

if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
  set(${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS C)
endif()

unset({CMAKE_FIND_PACKAGE_NAME}_REQUIRED_VARS)
foreach(COMP C CXX Fortran)
  if("${COMP}" IN_LIST ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS) 
    list(APPEND ${CMAKE_FIND_PACKAGE_NAME}_REQUIRED_VARS NETCDF_${COMP}_FOUND)
  endif()
endforeach()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(${CMAKE_FIND_PACKAGE_NAME} HANDLE_COMPONENTS REQUIRED_VARS ${CMAKE_FIND_PACKAGE_NAME}_REQUIRED_VARS)

