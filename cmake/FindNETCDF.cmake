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
	else()
		find_path(NETCDF_Fortran_INCLUDE_DIRECTORIES netcdf.inc HINTS $ENV{NETCDF_DIR}/include ENV NETCDF_Fortran_INCLUDE_DIRECTORIES)
		find_library(NETCDF_Fortran_LIBRARIES netcdff HINTS ${NETCDF_Fortran_INCLUDE_DIRECTORIES}/../lib)
	endif()	
endif()

# if one wants to link to the static NetCDF, also link these libraries for C: hdf5, hdf5_hl, curl
if(CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
	include(CheckFunctionExists)
	check_function_exists(nc__create HAVE_C_NETCDF)

	if(HAVE_C_NETCDF)
		set(NETCDF_C_INCLUDE_DIRECTORIES "")
		set(NETCDF_C_LIBRARIES "")
	else()
		find_path(NETCDF_C_INCLUDE_DIRECTORIES netcdf.h HINTS $ENV{NETCDF_DIR}/include ENV NETCDF_C_INCLUDE_DIRECTORIES)
		find_library(NETCDF_C_LIBRARIES netcdf HINTS ${NETCDF_C_INCLUDE_DIRECTORIES}/../lib)
	endif()
endif()

if(CMAKE_CXX_COMPILER_LOADED)
	find_path(NETCDF_CXX_INCLUDE_DIRECTORIES netcdf HINTS $ENV{NETCDF_DIR}/include ENV NETCDF_CXX_INCLUDE_DIRECTORIES)
	# the cray toolchain (e.g. hlrn) disables dynamic linking by default. to enable it at build time do e.g. "CRAYPE_LINK_TYPE=dynamic make".	
	find_library(NETCDF_CXX_LIBRARIES NAMES netcdf_c++4 netcdf-cxx4 HINTS ${NETCDF_CXX_INCLUDE_DIRECTORIES}/../lib)
	if(NETCDF_CXX_INCLUDE_DIRECTORIES AND NETCDF_C_INCLUDE_DIRECTORIES)
		list(APPEND NETCDF_CXX_INCLUDE_DIRECTORIES ${NETCDF_C_INCLUDE_DIRECTORIES})
	endif()
	if(NETCDF_CXX_LIBRARIES AND NETCDF_C_LIBRARIES)
		list(APPEND NETCDF_CXX_LIBRARIES ${NETCDF_C_LIBRARIES})
	endif()
endif()