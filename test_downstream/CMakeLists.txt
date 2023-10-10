cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(test_downstream LANGUAGES Fortran)

find_package(fesom REQUIRED)

message("fesom found: ${fesom_DIR}")
message("fesom_VERSION: ${fesom_VERSION}")
message("fesom_HAVE_OPENMP: ${fesom_HAVE_OPENMP}")
message("fesom_HAVE_IFS_INTERFACE: ${fesom_HAVE_IFS_INTERFACE}")

if(NOT TARGET fesom)
  message(FATAL_ERROR "fesom was installed wrongly: fesom target not found")
endif()

if (ASSERT_FESOM_DIR)
  if (NOT fesom_DIR MATCHES "${ASSERT_FESOM_DIR}")
    message(FATAL_ERROR "fesom was not found in expected location ${ASSERT_FESOM_DIR}")
  endif()
endif()

if (DEFINED ASSERT_HAVE_IFS_INTERFACE)
  if (ASSERT_HAVE_IFS_INTERFACE AND NOT fesom_HAVE_IFS_INTERFACE)
    message(FATAL_ERROR "fesom was installed wrongly: fesom_HAVE_IFS_INTERFACE=${fesom_HAVE_IFS_INTERFACE}")
  elseif (NOT ASSERT_HAVE_IFS_INTERFACE AND fesom_HAVE_IFS_INTERFACE)
    message(FATAL_ERROR "fesom was installed wrongly: fesom_HAVE_IFS_INTERFACE=${fesom_HAVE_IFS_INTERFACE}")
  endif()
endif()

get_target_property(target_type fesom TYPE)
if( NOT target_type MATCHES "STATIC_LIBRARY|SHARED_LIBRARY")
  message(FATAL_ERROR "fesom target's TYPE is expected to be STATIC_LIBRARY or SHARED_LIBRARY. Actual TYPE: ${target_type}")
endif ()

add_executable( main main.F90 )
target_link_libraries( main fesom )

