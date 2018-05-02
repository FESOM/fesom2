######################################################
# Makefile                
######################################################

######################################################
# Include architecture-specific definitions

include Makefile.in

###### Objects for Mesh Partitioning ################################################
# modules
MOD_INI =  fort_part.o  oce_modules.o  gen_modules_config.o gen_modules_partitioning.o gen_modules_rotate_grid.o  

OBJ_INI =  fvom_init.o \
           oce_local.o \
           gen_comm.o

# objects

MODULES = oce_modules.o \
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
          gen_modules_diag.o \
          psolve.o \
          oce_ale_mixing_pp.o \
          oce_tracer_mod.o \
	  fort_part.o \
          gen_input.o \
          io_meandata.o \
          io_restart.o \
          io_blowup.o \
          io_mesh_info.o

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
          oce_ale_fct_3d_adv.o 	\
          oce_ice_init_state.o 	\
          oce_shortwave_pene.o 	\
          oce_ale.o 		\
          oce_ale_tracer.o 	\
          ice_evp.o             \
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
          write_step_info.o



# oce_pressure_bv.o \
#          oce_redi_gm.o \
# Name of executables

EXE     = fesom.x
EXE_INI = fvom_ini.x

# TARGETS

default: run

run:  $(MODULES) $(OBJECTS)
	@echo "======= Building FESOM =========="
	$(LD) $(OPT) -o $(EXE) $(FOS_INC) $(MODULES) $(OBJECTS) \
	 $(MPI_LIB) $(LIB_LAP) $(LIB_PARMS) $(LIB_METIS) $(NCLIB) 
	cp -pf $(EXE) ../bin/.

run_ini: CPP_DEFS+=-DFVOM_INIT
run_ini: cleanomod $(MOD_INI) $(OBJ_INI)
	@echo "======= Building FESOM paritioning program =========="
	$(LD)  $(OPT) -o $(EXE_INI)  $(MOD_INI) $(OBJ_INI) \
	$(MPI_LIB) $(LIB_METIS) $(NCLIB)
	cp -pf $(EXE_INI) ../bin/.

.SUFFIXES: .c .F90 .o 

.c.o :
	$(CC) $(COPT) $(METIS_DEFS) $(CPP_DEFS) $(PARMS_DEFS) $(METIS_INC) $(PARMS_INC) $(MPI_INCLUDE) -c $*.c	

#.F90.o :
#	@echo $@
#	$(CPP) $(CPP_DEFS) $(CPP_SOL) $(PETSC_INC) $(PARMS_INC) $(PETSCCONF_INC) $(FOS_INC) $(MPI_INC) $(NCINC) $(LIBS_SLV) < $*.F90 > $*.f90
#	$(FC) $(OPT) $(CPP_SOL) $(CPP_DEFS) $(OASIS3_INC) $(FOS_INC)\
#	$(PARMS_INC) $(MPI_INCLUDE) $(NCINC)  -c $*.f90
.F90.o :
	@echo $@
	$(FC) $(CPP_DEFS) $(CPP_SOL) $(PETSC_INC) $(PARMS_INC) $(PETSCCONF_INC) $(FOS_INC) $(MPI_INC) $(NCINC) $(LIBS_SLV) $(OPT) $(CPP_SOL) $(CPP_DEFS) $(OASIS3_INC) $(FOS_INC)\
	$(PARMS_INC) $(MPI_INCLUDE) $(NCINC)  -c $*.F90

clean :
	rm -f *.o *.mod *~ *.f90 fvom.x fvom_ini.x

cleanomod:
	rm -f *.o *.mod
	CPP_DEFS=; export CPP_DEFS
