


# BLAS_<lang>_FOUND
# BLAS_<lang>_LIBRARIES
# BLAS_<lang>_COMPILER_WRAPPER


if(CMAKE_C_COMPILER_LOADED)
	#check_function_exists(cblas_sgemm SGEMM_C_FOUND) # this is also found if it is in the standard library paths, so we can not be sure we already have it in the wrapper

	include(CheckCSourceCompiles)
	set(CMAKE_REQUIRED_FLAGS -c)
	check_c_source_compiles([[
	#include <cblas.h>
	typedef void (*foo)(const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE,
		const enum CBLAS_TRANSPOSE, const int, const int,
		const int, const float, const float *,
		const int, const float *, const int,
		const float, float *, const int);
	foo func = &cblas_sgemm;
]] BLAS_C_COMPILER_WRAPPER)

	if(NOT BLAS_C_COMPILER_WRAPPER)
		# try to locate BLAS via the find module shipped with cmake
		# the module does not search for C based blas if the Fortran language is enabled
		include(${CMAKE_ROOT}/Modules/FindBLAS.cmake)
		set(BLAS_C_LIBRARIES ${BLAS_LIBRARIES})
		# unset the results from the original FindBlas.cmake
		unset(BLAS_FOUND)
		unset(BLAS_LIBRARIES)
	else()
		set(BLAS_C_LIBRARIES "") # we do not want to explicitly add libraries as the wrapper will do it
		set(BLAS_C_FOUND TRUE)
	endif()
	
	if(BLAS_C_LIBRARIES)
		set(BLAS_C_FOUND TRUE)
	endif()
endif()


# include(${CMAKE_ROOT}/Modules/FindBLAS.cmake)
# This module sets the following variables:
#
# ::
#
#   BLAS_FOUND - set to true if a library implementing the BLAS interface
#     is found
#   BLAS_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#     and -L).
#   BLAS_LIBRARIES - uncached list of libraries (using full path name) to
#     link against to use BLAS
#   BLAS95_LIBRARIES - uncached list of libraries (using full path name)
#     to link against to use BLAS95 interface
#   BLAS95_FOUND - set to true if a library implementing the BLAS f95 interface
#     is found
#   BLA_STATIC  if set on this determines what kind of linkage we do (static)
#   BLA_VENDOR  if set checks only the specified vendor, if not set checks
#      all the possibilities
#   BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
