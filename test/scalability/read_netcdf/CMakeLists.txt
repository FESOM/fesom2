cmake_minimum_required(VERSION 3.4)

project(read_netcdf Fortran)

# get our source files
file(GLOB sources_Fortran ${CMAKE_CURRENT_LIST_DIR}/*.F90)

include(${CMAKE_CURRENT_LIST_DIR}/../../../cmake/FindNETCDF.cmake)

add_executable(${PROJECT_NAME} ${sources_Fortran}
                               ${CMAKE_CURRENT_LIST_DIR}/../../../src/io_netcdf_module.F90
                               ${CMAKE_CURRENT_LIST_DIR}/../../../src/io_netcdf_nf_interface.F90
                               ${CMAKE_CURRENT_LIST_DIR}/../../../src/io_netcdf_file_module.F90
                               ${CMAKE_CURRENT_LIST_DIR}/../../../src/io_netcdf_attribute_module.F90
                               )

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero -no-wrap-margin)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	target_compile_options(${PROJECT_NAME} PRIVATE -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none)
	if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
		target_compile_options(${PROJECT_NAME} PRIVATE -fallow-argument-mismatch) # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
	endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray )
       target_compile_options(${PROJECT_NAME} PRIVATE -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64)
endif()

target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
