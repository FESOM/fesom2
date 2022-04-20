# Todo: set the following flags:
#CCFLAGS = -Wall -Wextra -lstdc++ -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib
#LDFLAGS = -lstdc++ -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib

cmake_minimum_required(VERSION 3.4)

message("Generating hecuba interface library" )
#  the FortranCInterface requires the C language to be enabled
project(io_hecuba_cpp CXX C Fortran)

if(DOCKER_HECUBA)
    message("Setting path to Hecuba includes for Docker container platform" )
    set( HECUBA_ROOT "/usr/local/" )
else()
    message("Setting path to Hecuba includes for Marenostrum IV platform")
    set( HECUBA_ROOT "/apps/HECUBA/1.1" )    
endif()

# set Hecuba compilation flags - c++
set(HECUBA_CCFLAGS -Wall -Wextra -std=c++11 -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib)

# set Hecuba compilation flags - c++
# set(HECUBA_INCLUDES -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib)

# get our source files
file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(HecubaWrapperFortranMacros.h MACRO_NAMESPACE "HecubaWrapperFortranMacros_" SYMBOLS init_datamodel_c hecuba_output_c)

add_library(${PROJECT_NAME} ${sources_CXX})

target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_LIST_DIR}
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)

target_compile_options(${PROJECT_NAME} PRIVATE ${HECUBA_CCFLAGS})
