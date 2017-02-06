.SUFFIXES: .F90 .f90 .c .o

# Model Object Files
OBJECTS = oce_read_mesh.o \
	interp.o 	\
	interp_el.o 	\
	make_grid_info.o \
	make_diag_uv_norm3.o \
	make_diag_moc_w.o \
	make_diag_uv_curl3.o \
	make_diag_ts3.o \
	do_work.o

MODULES = oce_modules.o             \
	  gen_modules_config.o      \
	  gen_modules_rotate_grid.o \
	  mesh_setup.o              \
	  g_oce_2_reg.o             \
          oce_read_data.o

# NETCDF
NCINC = $(NETCDF_DIR)/include
NCLIB = $(NETCDF_DIR)/lib

EXE = do_work.x

CC = cc
CCFLAGS = -g
FC = ifort #/opt/intel/fc/9.1.033/bin/ifort
FC_OPT = -r8 -i4 -g -fpp -bc -traceback
LD = $(FC)

CPP = cpp

$(EXE): $(MODULES) $(OBJECTS)
	$(LD) -o $(EXE) -L$(NCLIB) $(MODULES) $(OBJECTS) -lnetcdff -lnetcdf

.F90.o :
	$(CPP) -P -C $(CPP_DEFS) -I$(NCINC) < $*.F90 > $*.f90
	$(FC) $(FC_OPT) -c $*.f90 #-I$(INC_SLV) -I$(NCINC)
	rm -f  $*.f90

.f.o :
	$(FC) -c $*.f

.c.o:
	$(CC) -c $(CCFLAGS) $*.c

clean :
	rm -f *.o *.mod *~ *.f90 *.x

run: $(EXE)
	./$(EXE)
