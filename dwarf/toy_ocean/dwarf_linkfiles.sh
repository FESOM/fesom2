#!/bin/bash
#_______________________________________________________________________________
# link environment variables and configure.sh files from the main repository level
ln -sf ../../env env
ln -sf ../../env.sh env.sh
ln -sf ../../configure.sh configure.sh

#_______________________________________________________________________________
# Create bin directory for env.sh (needed for current_shell_path)
mkdir -p bin

#_______________________________________________________________________________
# create local source folder for the dwarf
if [ -d  "src/" ] ; then rm -r src/ ; fi
mkdir src/
cd src/

#_______________________________________________________________________________
# link main dwarf files to the local src/ folder
ln -sf ../dwarf_ini/fesom_toy_main.F90 fesom_toy_main.F90
ln -sf ../dwarf_ini/fesom_toy_module.F90 fesom_toy_module.F90
ln -sf ../dwarf_ini/CMakeLists.txt CMakeLists.txt
ln -sf ../dwarf_ini/fesom_toy-config.cmake.in fesom_toy-config.cmake.in

#_______________________________________________________________________________
export which_path="../../../src/"

# for downloading from specific github branch replace ln -s with wget
# export which_branch=main
# export which_path=https://raw.githubusercontent.com/FESOM/fesom2/${which_branch}/src

#_______________________________________________________________________________
# Core Modules (MOD_*)
export which_files="MOD_DYN.F90
                    MOD_ICE.F90
                    MOD_MESH.F90
                    MOD_PARTIT.F90
                    MOD_READ_BINARY_ARRAYS.F90
                    MOD_TRACER.F90
                    MOD_WRITE_BINARY_ARRAYS.F90"

# Ocean Modules (oce_*)
export which_files="${which_files}
                    oce_adv_tra_driver.F90
                    oce_adv_tra_fct.F90
                    oce_adv_tra_hor.F90
                    oce_adv_tra_ver.F90
                    oce_ale.F90
                    oce_ale_mixing_kpp.F90
                    oce_ale_mixing_pp.F90
                    oce_ale_pressure_bv.F90
                    oce_ale_ssh_splitexpl_subcycl.F90
                    oce_ale_tracer.F90
                    oce_ale_vel_rhs.F90
                    oce_dyn.F90
                    oce_fer_gm.F90
                    oce_mesh.F90
                    oce_mo_conv.F90
                    oce_modules.F90
                    oce_muscl_adv.F90
                    oce_setup_step.F90
                    oce_shortwave_pene_dbgyre.F90
                    oce_spp.F90
                    oce_tracer_mod.F90"

# Global Modules (gen_* and cavity_*)
export which_files="${which_files}
                    cavity_param.F90
                    gen_events.F90
                    gen_forcing_init.F90
                    gen_halo_exchange.F90
                    gen_ic3d.F90
                    gen_interpolation.F90
                    gen_model_setup.F90
                    gen_modules_backscatter.F90
                    gen_modules_clock.F90
                    gen_modules_config.F90
                    gen_modules_diag.F90
                    gen_modules_forcing.F90
                    gen_modules_partitioning.F90
                    gen_modules_read_NetCDF.F90
                    gen_modules_rotate_grid.F90
                    gen_support.F90"

# I/O Modules (io_*)
export which_files="${which_files}
                    io_blowup.F90
                    io_data_strategy.F90
                    io_fesom_file.F90
                    io_gather.F90
                    io_meandata.F90
                    io_mesh_info.F90
                    io_netcdf_attribute_module.F90
                    io_netcdf_file_module.F90
                    io_netcdf_workaround_module.F90
                    io_restart.F90
                    io_restart_derivedtype.F90
                    io_restart_file_group.F90
                    io_scatter.F90"

# Forcing Modules
export which_files="${which_files}
                    forcing_lookahead_reader_module.F90
                    forcing_provider_async_module.F90
                    forcing_provider_netcdf_module.F90
                    gen_surface_forcing.F90"

# Toy-Specific Modules
export which_files="${which_files}
                    toy_channel_dbgyre.F90
                    toy_channel_soufflet.F90"

# Iceberg Modules (needed for compilation even though not actively used)
export which_files="${which_files}
                    icb_allocate.F90
                    icb_coupling.F90
                    icb_dyn.F90
                    icb_elem.F90
                    icb_modules.F90
                    icb_step.F90
                    icb_thermo.F90"

# Support Modules
export which_files="${which_files}
                    async_threads_module.F90
                    command_line_options.F90
                    fortran_utils.F90
                    info_module.F90
                    mod_transit.F90
                    mpi_topology_module.F90
                    nvfortran_subarray_workaround.F90
                    solver.F90
                    write_step_info.F90"

# CVMix drivers (conditional, will be linked if CVMIX is enabled)
export which_files="${which_files}
                    cvmix_driver/cvmix_idemix.F90
                    cvmix_driver/cvmix_kinds_and_types_addon.F90
                    cvmix_driver/cvmix_tke.F90
                    cvmix_driver/cvmix_utils_addon.F90
                    cvmix_driver/gen_modules_cvmix_idemix.F90
                    cvmix_driver/gen_modules_cvmix_kpp.F90
                    cvmix_driver/gen_modules_cvmix_pp.F90
                    cvmix_driver/gen_modules_cvmix_tidal.F90
                    cvmix_driver/gen_modules_cvmix_tke.F90"

# IFS interface (for MULTIO support)
export which_files="${which_files}
                    ifs_interface/iom.F90
                    ifs_interface/mpp_io.F90"

# Version info
export which_files="${which_files}
                    fesom_version_info.F90"

# Profiler (optional)
export which_files="${which_files}
                    fesom_profiler.F90"

#_______________________________________________________________________________
# link the necessary main src files to local src directory
for file in ${which_files}; do
    # Create subdirectories if needed (for cvmix_driver, ifs_interface)
    if [[ $file == *"/"* ]]; then
        mkdir -p $(dirname ${file})
    fi
    ln -sf ${which_path}/${file} ${file}
    # wget ${which_path}/${file}
    # cp ${which_path}/${file} .
done

#_______________________________________________________________________________
# link async_threads_cpp directory for C++ threading support
ln -sf ${which_path}/async_threads_cpp async_threads_cpp

#_______________________________________________________________________________
cd ../

echo "Linkfiles created successfully!"
echo "You can now build the toy ocean model by running:"
echo "  ./configure.sh <platform> [options]"
