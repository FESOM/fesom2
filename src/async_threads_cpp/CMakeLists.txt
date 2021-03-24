cmake_minimum_required(VERSION 3.4)

project(async_threads_cpp CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(ThreadsManagerFCMacros.h MACRO_NAMESPACE "ThreadsManagerFCMacros_" SYMBOLS init_ccall begin_ccall end_ccall)

add_library(${PROJECT_NAME} ${sources_CXX})
target_include_directories(${PROJECT_NAME}
	INTERFACE ${CMAKE_CURRENT_LIST_DIR}
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package( Threads )
target_link_libraries( ${PROJECT_NAME} Threads::Threads)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 11)

#check_symbol_exists(<symbol> <files> <variable>)
include(CheckSymbolExists)
# https://bugs.llvm.org/show_bug.cgi?id=12730
check_symbol_exists(__GCC_HAVE_SYNC_COMPARE_AND_SWAP_1 "stdlib.h" HAVE_SWAP_1)
check_symbol_exists(__GCC_HAVE_SYNC_COMPARE_AND_SWAP_2 "stdlib.h" HAVE_SWAP_2)
check_symbol_exists(__GCC_HAVE_SYNC_COMPARE_AND_SWAP_4 "stdlib.h" HAVE_SWAP_4)
check_symbol_exists(__GCC_HAVE_SYNC_COMPARE_AND_SWAP_8 "stdlib.h" HAVE_SWAP_8)

target_compile_options(${PROJECT_NAME} PRIVATE
	$<$<NOT:$<BOOL:${HAVE_SWAP_1}>>:-D__GCC_HAVE_SYNC_COMPARE_AND_SWAP_1>
	$<$<NOT:$<BOOL:${HAVE_SWAP_2}>>:-D__GCC_HAVE_SYNC_COMPARE_AND_SWAP_2>
	$<$<NOT:$<BOOL:${HAVE_SWAP_4}>>:-D__GCC_HAVE_SYNC_COMPARE_AND_SWAP_4>
	$<$<NOT:$<BOOL:${HAVE_SWAP_8}>>:-D__GCC_HAVE_SYNC_COMPARE_AND_SWAP_8>
)
