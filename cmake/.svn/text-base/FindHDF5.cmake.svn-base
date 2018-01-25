# HDF5_C_INCLUDE_DIRECTORIES
# HDF5_C_LIBRARIES


if(CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
   include(CheckFunctionExists)
   check_function_exists(H5get_libversion HAVE_C_HDF5)
   if(HAVE_C_HDF5)
      set(HDF5_C_INCLUDE_DIRECTORIES "")
      set(HDF5_C_LIBRARIES "")
   else()
      find_path(HDF5_C_INCLUDE_DIRECTORIES hdf5.h HINTS $ENV{HDF5_DIR}/include ENV HDF5_C_INCLUDE_DIRECTORIES)
      find_library(HDF5_C_LIBRARIES hdf5 HINTS ${HDF5_C_INCLUDE_DIRECTORIES}/../lib)
   endif()
endif()
