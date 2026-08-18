# NETCDF_Fortran_INCLUDE_DIRECTORIES
# NETCDF_Fortran_LIBRARIES
# NETCDF_C_INCLUDE_DIRECTORIES
# NETCDF_C_LIBRARIES

# Enhanced FindNETCDF with compiler compatibility checking
# If BUILD_NETCDF option is enabled and system NetCDF is incompatible,
# this will be skipped as the external build will set the variables directly

# Skip system NetCDF search if BUILD_NETCDF is enabled
if(BUILD_NETCDF)
    # Variables will be set by BuildNetCDF.cmake
    return()
endif()

# Function to test NetCDF compatibility with current compiler
function(test_netcdf_compatibility NETCDF_INCLUDE_DIR NETCDF_LIB COMPILER_WORKS_VAR)
    if(NOT NETCDF_INCLUDE_DIR OR NOT NETCDF_LIB)
        set(${COMPILER_WORKS_VAR} FALSE PARENT_SCOPE)
        return()
    endif()
    
    # Create a simple test program
    set(TEST_SOURCE_DIR ${CMAKE_BINARY_DIR}/netcdf_test)
    file(MAKE_DIRECTORY ${TEST_SOURCE_DIR})
    
    # Simple Fortran test
    file(WRITE ${TEST_SOURCE_DIR}/test_netcdf.F90
"program test_netcdf
    use netcdf
    implicit none
    integer :: ncid, status
    status = nf90_create('test.nc', NF90_CLOBBER, ncid)
    if (status /= NF90_NOERR) then
        stop 1
    endif
    status = nf90_close(ncid)
end program test_netcdf")
    
    # Try to compile and link the test
    try_compile(COMPILE_RESULT
        ${TEST_SOURCE_DIR}
        ${TEST_SOURCE_DIR}/test_netcdf.F90
        CMAKE_FLAGS
            "-DINCLUDE_DIRECTORIES:STRING=${NETCDF_INCLUDE_DIR}"
            "-DLINK_LIBRARIES:STRING=${NETCDF_LIB}"
        OUTPUT_VARIABLE COMPILE_OUTPUT
    )
    
    set(${COMPILER_WORKS_VAR} ${COMPILE_RESULT} PARENT_SCOPE)
    
    if(NOT COMPILE_RESULT)
        message(STATUS "NetCDF compatibility test failed with current compiler")
        message(STATUS "Compile output: ${COMPILE_OUTPUT}")
    endif()
endfunction()

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
			# Test NetCDF compatibility with current compiler
			test_netcdf_compatibility("${NETCDF_Fortran_INCLUDE_DIRECTORIES}" "${NETCDF_Fortran_LIBRARIES}" NETCDF_FORTRAN_COMPATIBLE)
			
			if(NETCDF_FORTRAN_COMPATIBLE)
				set(NETCDF_Fortran_FOUND 1)
				message(STATUS "System NetCDF-Fortran is compatible with ${CMAKE_Fortran_COMPILER_ID} compiler")
			else()
				message(WARNING "System NetCDF-Fortran found but incompatible with ${CMAKE_Fortran_COMPILER_ID} compiler")
				message(STATUS "Consider using -DBUILD_NETCDF=ON to build compatible NetCDF libraries")
				set(NETCDF_Fortran_FOUND 0)
			endif()
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

