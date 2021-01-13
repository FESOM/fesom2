cmake_minimum_required(VERSION 3.13) # 20-11-2018
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")

project(fesom C Fortran)

option(DISABLE_MULTITHREADING "disable asynchronous operations" OFF)

# generate a custom file from fesom_version_info.F90 which includes the current git SHA
set(FESOM_ORIGINAL_VERSION_FILE ${CMAKE_CURRENT_LIST_DIR}/fesom_version_info.F90)
set(FESOM_GENERATED_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/fesom_version_info-generated.F90)
add_custom_command(OUTPUT 5303B6F4_E4F4_45B2_A6E5_8E2B9FB5CDC4 ${FESOM_GENERATED_VERSION_FILE} # the first arg to OUTPUT is a name for a file we never create to make sure this command will run on every re-build (let our file be the second arg, as the first file is inadvertently removed by make)
                   COMMAND ${CMAKE_COMMAND} -DFESOM_ORIGINAL_VERSION_FILE=${FESOM_ORIGINAL_VERSION_FILE} -DFESOM_GENERATED_VERSION_FILE=${FESOM_GENERATED_VERSION_FILE} -P GitRepositoryInfo.cmake
                   WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                   COMMENT "determining fesom git SHA ...")

# depends on the metis library
#add_subdirectory(../lib/metis-5.1.0 ${PROJECT_BINARY_DIR}/metis)
#include_directories(../lib/metis-5.1.0/include)
# depends on the parms library
add_subdirectory(../lib/parms ${PROJECT_BINARY_DIR}/parms)

add_subdirectory(async_threads_cpp)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)

# create our binary (set its name to name of this project)
add_executable(fesom
  async_threads_module.F90
  cavity_param.F90
  cpl_driver.F90
  cvmix_background.F90
  cvmix_convection.F90
  cvmix_ddiff.F90
  cvmix_idemix.F90
  cvmix_kinds_and_types.F90
  cvmix_kpp.F90
  cvmix_math.F90
  cvmix_put_get.F90
  cvmix_shear.F90
  cvmix_tidal.F90
  cvmix_tke.F90
  cvmix_utils.F90
  fesom_version_info-generated.F90
  forcing_lookahead_reader_module.F90
  forcing_provider_async_module.F90
  forcing_provider_netcdf_module.F90
  fvom_main.F90
  gen_bulk_formulae.F90
  gen_comm.F90
  gen_events.F90
  gen_forcing_couple.F90
  gen_forcing_init.F90
  gen_halo_exchange.F90
  gen_ic3d.F90
  gen_interpolation.F90
  gen_model_setup.F90
  gen_modules_clock.F90
  gen_modules_config.F90
  gen_modules_cvmix_idemix.F90
  gen_modules_cvmix_kpp.F90
  gen_modules_cvmix_pp.F90
  gen_modules_cvmix_tidal.F90
  gen_modules_cvmix_tke.F90
  gen_modules_diag.F90
  gen_modules_forcing.F90
  gen_modules_gpot.F90
  gen_modules_partitioning.F90
  gen_modules_read_NetCDF.F90
  gen_modules_rotate_grid.F90
  gen_support.F90
  gen_surface_forcing.F90
  ice_EVP.F90
  ice_fct.F90
  ice_maEVP.F90
  ice_modules.F90
  ice_oce_coupling.F90
  ice_setup_step.F90
  ice_thermo_cpl.F90
  ice_thermo_oce.F90
  io_blowup.F90
  io_data_strategy.F90
  io_gather.F90
  io_meandata.F90
  io_mesh_info.F90
  io_netcdf_workaround_module.F90
  io_restart.F90
  MOD_MESH.F90
  mpi_topology_module.F90
  oce_adv_tra_driver.F90
  oce_adv_tra_fct.F90
  oce_adv_tra_hor.F90
  oce_adv_tra_ver.F90
  oce_ale.F90
  oce_ale_mixing_kpp.F90
  oce_ale_mixing_pp.F90
  oce_ale_pressure_bv.F90
  oce_ale_tracer.F90
  oce_ale_vel_rhs.F90
  oce_dyn.F90
  oce_fer_gm.F90
  oce_ice_init_state.F90
  oce_local.F90
  oce_mesh.F90
  oce_mo_conv.F90
  oce_modules.F90
  oce_muscl_adv.F90
  oce_setup_step.F90
  oce_shortwave_pene.F90
  oce_spp.F90
  oce_tracer_mod.F90
  oce_vel_rhs_vinv.F90
  psolve.c
  psolve_feom.c
  toy_channel_soufflet.F90
  write_step_info.F90
)
target_compile_definitions(fesom PRIVATE
  PARMS
  -DMETIS_VERSION=5 -DPART_WEIGHTED -DMETISRANDOMSEED=35243
  $<$<BOOL:${DISABLE_MULTITHREADING}>:DISABLE_MULTITHREADING>

  $<$<COMPILE_LANGUAGE:C>:
    USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2
   >
)

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
# https://cmake.org/cmake/help/latest/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_%3CLANG%3E_COMPILER_ID
target_compile_options(fesom PRIVATE
  $<$<BOOL:${FESOM_COUPLED}>:__oasis>
  $<$<BOOL:${OIFS_COUPLED}>:__oifs>
  $<$<BOOL:${VERBOSE}>:VERBOSE>

  # or use $<COMPILE_LANG_AND_ID:language,compiler_ids>
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:
    -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip -init=zero -no-wrap-margin
   >

  $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:
    -O3 -finit-local-zero  -finline-functions -march=native -fimplicit-none  -fdefault-real-8 -ffree-line-length-none
    $<$<VERSION_GREATER_EQUAL:${CMAKE_Fortran_COMPILER_VERSION},10>:
      -fallow-argument-mismatch # gfortran v10 is strict about erroneous API calls: "Rank mismatch between actual argument at (1) and actual argument at (2) (scalar and rank-1)"
     >
   >

  $<$<COMPILE_LANG_AND_ID:Fortran,Cray>:
    -c -emf -hbyteswapio -hflex_mp=conservative -hfp1 -hadd_paren -Ounroll0 -hipa0 -r am -s real64
   >

  $<$<COMPILE_LANG_AND_ID:C,PGI>:
    -fast -fastsse -O3 -Minline -Munroll
   >
  $<$<COMPILE_LANG_AND_ID:Fortran,PGI>:
   -fast -fastsse -O3 -Mallocatable=95 -Minline -Munroll -Mr8 -pgf90libs
    $<$<BOOL:${FESOM_OPENACC}>:
      -acc -ta=tesla:cc35
      -Minfo=accel
     >
   >
  )

target_link_options(fesom PRIVATE
  $<$<Fortran_COMPILER_ID:PGI>:
    -acc -ta=tesla:cc35
  >
)
if(${FESOM_COUPLED})
   if(${CRAY})
      include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS_cray.cmake)
   else()
      if(${OIFS_COUPLED})
        include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS_oifs.cmake)
      else()
        include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
      endif()
   endif()
endif()
target_include_directories(fesom PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(fesom PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(fesom PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(fesom ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES})
target_link_libraries(fesom ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
target_link_libraries(fesom parms) #metis
target_link_libraries(fesom async_threads_cpp)
set_target_properties(fesom PROPERTIES LINKER_LANGUAGE Fortran)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom.x" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/fesom DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
