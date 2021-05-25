######################################################
# Makefile                
######################################################

######################################################
# Include architecture-specific definitions

include ../Makefile.in

# Definition of metis include and library
METIS_DIR = ../lib/metis-5.1.0
LIB_METIS = -L$(METIS_DIR)/lib -lmetis
METIS_INC = -I$(METIS_DIR)/include

# Definition of pARMS include and library
PARMS_DIR = ../lib/parms
PARMS_INC = -I$(PARMS_DIR)/include -I$(PARMS_DIR)/src/include
LIB_PARMS = -L$(PARMS_DIR)/lib -lparms
CPP_SOL   = -DPARMS

###### Objects for Mesh Partitioning ################################################
# modules
MOD_INI =  fort_part.o  oce_modules.o MOD_MESH.o gen_modules_config.o gen_modules_partitioning.o gen_modules_rotate_grid.o  

OBJ_INI =  fvom_init.o \
           oce_local.o \
           gen_comm.o

# objects

MODULES = oce_modules.o \
          info_module.o \
          command_line_options.o \
          MOD_MESH.o \
          ice_modules.o \
          gen_modules_config.o \
          gen_modules_partitioning.o \
          gen_modules_clock.o \
          gen_modules_rotate_grid.o \
          gen_modules_read_NetCDF.o \
          gen_modules_forcing.o \
          gen_halo_exchange.o \
          gen_support.o \
          oce_ale_mixing_kpp.o \
          oce_adv_tra_hor.o \
          oce_adv_tra_ver.o \
          oce_adv_tra_fct.o \
          oce_adv_tra_driver.o \
          gen_modules_diag.o \
          psolve.o \
          oce_ale_mixing_pp.o \
          oce_tracer_mod.o \
          cvmix_kinds_and_types.o \
          cvmix_utils.o \
          cvmix_put_get.o \
          cvmix_tke.o \
	  cvmix_idemix.o \
          gen_modules_cvmix_idemix.o \
          gen_modules_cvmix_tke.o \
          cvmix_math.o \
          cvmix_kpp.o \
          gen_modules_cvmix_kpp.o \
          cvmix_tidal.o \
          gen_modules_cvmix_tidal.o \
          cvmix_shear.o \
          gen_modules_cvmix_pp.o \
          async_threads_module.o \
          forcing_provider_netcdf_module.o \
          forcing_lookahead_reader_module.o \
          forcing_provider_async_module.o \
          io_gather.o \
          mpi_topology_module.o \
          io_netcdf_workaround_module.o \
          io_data_strategy.o \
          fesom_version_info.o \
          io_meandata.o \
          io_restart.o \
          io_blowup.o \
          io_mesh_info.o \
          gen_ic3d.o \
          gen_surface_forcing.o \
          gen_modules_gpot.o \
          toy_channel_soufflet.o 

OBJECTS=  fvom_main.o 		\
          gen_comm.o 		\
          oce_setup_step.o 	\
          oce_mesh.o 		\
          oce_dyn.o 		\
          oce_ale_vel_rhs.o 	\
          oce_vel_rhs_vinv.o 	\
	  oce_ale_pressure_bv.o \
	  oce_fer_gm.o 		\
          oce_muscl_adv.o 	\
          oce_ice_init_state.o 	\
          oce_shortwave_pene.o 	\
          oce_ale.o 		\
          oce_ale_tracer.o 	\
          cavity_param.o        \
          ice_EVP.o             \
          ice_maEVP.o 		\
          ice_setup_step.o 	\
          ice_fct.o 		\
          ice_oce_coupling.o 	\
          ice_thermo_oce.o      \
          gen_model_setup.o 	\
          gen_forcing_init.o 	\
          gen_bulk_formulae.o 	\
          gen_forcing_couple.o 	\
          gen_interpolation.o 	\
          gen_events.o 		\
          write_step_info.o     \
          oce_mo_conv.o         \
          oce_spp.o  


# oce_pressure_bv.o \
#          oce_redi_gm.o \
# Name of executables

EXE     = fesom.x
EXE_INI = fesom_ini.x

# TARGETS

default: run

run:  $(MODULES) $(OBJECTS)
	@echo "======= Building FESOM =========="
	$(LD) $(OPT) -o $(EXE) $(MODULES) $(OBJECTS) \
	 $(MPI_LIB) $(LIB_LAP) $(LIB_PARMS) $(NC_LIB) 
#	cp -pf $(EXE) ../bin/.

run_ini: CPP_DEFS+=-DFVOM_INIT
run_ini: cleanomod $(MOD_INI) $(OBJ_INI)
	@echo "======= Building FESOM partioning program =========="
	$(LD)  $(OPT) -o $(EXE_INI)  $(MOD_INI) $(OBJ_INI) \
	$(MPI_LIB) $(LIB_METIS) $(NC_LIB)
	cp -pf $(EXE_INI) ../bin/.

.SUFFIXES: .c .F90 .o 

.c.o :
	$(MPICC) $(COPT) $(CPP_SOL) $(METIS_DEFS) $(PARMS_DEFS) $(CPP_DEFS) $(METIS_INC) $(PARMS_INC) $(MPI_INC) -c $*.c	

.F90.o :
	@echo $@
	$(MPIFC) $(CPP_DEFS) $(CPP_SOL) $(NC_INC) $(OPT) $(OASIS3_INC) $(PARMS_INC) $(MPI_INC)  -c $*.F90

clean :
	rm -f *.o *.mod *~ *.f90 $(EXE) $(EXE_INI)

cleanomod:
	rm -f *.o *.mod
	CPP_DEFS=; export CPP_DEFS
