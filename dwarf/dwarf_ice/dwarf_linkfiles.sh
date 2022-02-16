#!/bin/bash
#_______________________________________________________________________________
# link environment variables and configure.sh files from the main repository level
ln -s ../../env env
ln -s ../../env.sh env.sh
ln -s ../../configure.sh configure.sh

#_______________________________________________________________________________
# create local source folder for the dwarf
if [ -d  "src/" ] ; then rm -r src/ ; fi 
mkdir src/
cd src/

#_______________________________________________________________________________
# link main dwarf files to the local src/ folder
ln -s ../dwarf_ini/fesom.F90 fesom.F90
ln -s ../dwarf_ini/CMakeLists.txt CMakeLists.txt

#_______________________________________________________________________________
export which_path="../../../src/"

# for downloading from specific github branch replace ln -s with wget
# export which_branch=refactoring
# export which_branch=refactoring_dwarf_ice
# export which_path=https://raw.githubusercontent.com/FESOM/fesom2/${which_branch}/src

#_______________________________________________________________________________
export which_files="associate_mesh_def.h 
                    associate_mesh_ass.h 
                    associate_part_def.h 
                    associate_part_ass.h
                    MOD_MESH.F90
                    MOD_PARTIT.F90  
                    MOD_TRACER.F90
                    MOD_DYN.F90
                    MOD_ICE.F90
                    MOD_READ_BINARY_ARRAYS.F90
                    MOD_WRITE_BINARY_ARRAYS.F90
                    ice_modules.F90
                    ice_EVP.F90
                    ice_maEVP.F90
                    ice_fct.F90
                    gen_halo_exchange.F90
                    gen_modules_partitioning.F90
                    gen_modules_config.F90
                    io_restart_derivedtype.F90
                    fortran_utils.F90
                    oce_modules.F90
                    "
#_______________________________________________________________________________
# link the ther necessary main src files to local src directory 
for file in ${which_files}; do
    ln -s ${which_path}/${file} ${file}
    # wget ${which_path}/${file}
    # cp ${which_path}/${file} .
done

#_______________________________________________________________________________
cd ../

