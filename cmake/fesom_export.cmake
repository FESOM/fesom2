
macro( fesom_export )

set( options )
set( single_value_args )
set( multi_value_args  TARGETS )

cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

set( PROJECT_TARGETS_FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake" )
file( REMOVE ${PROJECT_TARGETS_FILE} )


foreach( tgt ${_PAR_TARGETS} )
   install( TARGETS ${tgt}
      EXPORT  ${PROJECT_NAME}-targets
      RUNTIME DESTINATION ${INSTALL_BIN_DIR}
      LIBRARY DESTINATION ${INSTALL_LIB_DIR}
      ARCHIVE DESTINATION ${INSTALL_LIB_DIR} )
   set_target_properties( ${tgt} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
   export( TARGETS ${tgt} APPEND FILE "${PROJECT_TARGETS_FILE}" )
endforeach()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_NAME}-config-version.cmake"
    VERSION ${${PROJECT_NAME}_VERSION}
    COMPATIBILITY AnyNewerVersion)

if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake.in )
   configure_package_config_file(${PROJECT_NAME}-config.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
      INSTALL_DESTINATION ${INSTALL_CMAKE_DIR})
else()
   file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake.in "include(${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}-targets.cmake)")
   configure_package_config_file(${PROJECT_NAME}-config.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
      INSTALL_DESTINATION ${INSTALL_CMAKE_DIR})
endif()

install(EXPORT ${PROJECT_NAME}-targets DESTINATION "${INSTALL_CMAKE_DIR}")
install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake
        DESTINATION "${INSTALL_CMAKE_DIR}")

install( DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR}
         DESTINATION module/${PROJECT_NAME}
         COMPONENT modules )

# Define ${PROJECT_NAME}_DIR in PARENT_SCOPE so that a `find_package( <this-project> )` in a bundle
# will easily find the project without requiring a `HINT <this-project>_BINARY_DIR` argument [ECBUILD-460]
if( NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR )
    # Guard needed because PARENT_SCOPE cannot be used in top-level CMake project

    set( ${PROJECT_NAME}_DIR ${PROJECT_BINARY_DIR} PARENT_SCOPE )
endif()

endmacro()
