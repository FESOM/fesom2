cmake_minimum_required(VERSION 3.4)

project(io_hecuba CXX) # currently not generating fortran lib automatically

# get our source files
file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

add_library(${PROJECT_NAME} SHARED ${sources_CXX})

set(HECUBA_ROOT $ENV{HECUBA_ROOT})
# TODO: check if HECUBA_ROOT exists if not atleast raise

message("HECUBA ROOT is ${HECUBA_ROOT}")
target_include_directories(${PROJECT_NAME}
        INTERFACE ${CMAKE_CURRENT_LIST_DIR}
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
	PUBLIC "${HECUBA_ROOT}/include"
)

#using public below works but will libs in $HECUBA_ROOT may conflict with others
set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-Wl,-rpath,${HECUBA_ROOT}/lib/")
#set_target_properties(${PROJECT_NAME} PROPERTIES CMAKE_INSTALL_RPATH_USE_LINK_PATH "${HECUBA_ROOT}/lib")

#target_link_directories(${PROJECT_NAME} PUBLIC "${HECUBA_ROOT}/lib")
#target_link_libraries(${PROJECT_NAME} PUBLIC hfetch)

target_link_libraries(${PROJECT_NAME} PUBLIC -L${HECUBA_ROOT}/lib -lhfetch)

#target_link_libraries(${PROJECT_NAME}_C -L${LIBDILL_LIB} -ldill -pthread -lssl -lcrypto)


#target_compile_options(${PROJECT_NAME} PRIVATE -std=c++11 -fpic -shared)
