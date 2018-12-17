# normally this files  need not be edited -- architecture-dependent
# make commands are in makefile.in -- 
# this makefile will make the library -- for testing PARARMS -- go to
# examples/grid or examples/general and see README files..

include ../../Makefile.in

# path of the header files of pARMS
ICFLAGS			=	-I./include $(MPI_INC)

# path of the header files for implmentaion of pARMS
ISRCINC	                = -I./src/include	

# path of the header files
IFFLAGS			=  -I./include

# library path and directory declaration
LIB	=    ./lib/libparms.a
DIRS    =    ./ ./include ./lib ./src ./src/include ./src/DDPQ 

# files from the src directory
OBJ1	=    ./src/parms_comm.o  \
	     ./src/parms_map.o  ./src/parms_mat.o  \
	     ./src/parms_mat_vcsr.o ./src/parms_mat_dvcsr.o  \
             ./src/parms_mem.o ./src/parms_operator.o  \
             ./src/parms_ilu_vcsr.o \
	     ./src/parms_pc.o \
	     ./src/parms_pc_bj.o  ./src/parms_pc_ras.o  \
	     ./src/parms_pc_schur.o ./src/parms_pc_schurras.o \
             ./src/parms_qsplit.o ./src/parms_solver.o  \
             ./src/parms_table.o ./src/parms_timer.o  \
             ./src/parms_vec.o ./src/parms_viewer.o \
	     ./src/fgmres.o ./src/bicgstab.o ./src/cg.o ./src/pbicgstab.o \
	     ./src/pbicgstab_ras.o ./src/bicgstab_ras.o\
             ./src/gmres.o ./src/parms_complex.o 
#$(HYPRE_OBJ) \

# files from DDPQ directory
OBJ2 	=    ./src/DDPQ/arms2.o ./src/DDPQ/ilutpC.o  \
	     ./src/DDPQ/MatOps.o ./src/DDPQ/misc.o  \
	     ./src/DDPQ/PQ.o ./src/DDPQ/piluNEW.o  \
             ./src/DDPQ/setblks.o ./src/DDPQ/sets.o  \
	     ./src/DDPQ/svdInvC.o 

OBJ     =     $(OBJ1) $(OBJ2)

default: $(LIB) 

all: $(LIB) tests

$(LIB): $(OBJ)
	if [ ! -d lib ]; then \
	  mkdir lib; \
	fi
	$(AR) $(LIB) $(OBJ) 

.c.o:
	${MPICC} ${ICFLAGS} ${ISRCINC} $(COPT) \
	${PARMS_DEFS} $< -c -o $@

cleanall: 
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir;  rm -rf *.a *.o *.ex* *core* out* \#*  *~) ;\
          done

cleanobj:
	@for dir in $(DIRS) ;\
	do \
          echo cleaning $$dir ;\
          (cd $$dir;  rm -rf  *.o  *core* out* sol.* \#*  *~) ;\
          done

