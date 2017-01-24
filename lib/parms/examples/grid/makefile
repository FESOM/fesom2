SHELL	= /bin/sh

include ../../makefile.in

# path of the header files of pARMS
ICFLAGS			=	-I../../include

# path of the header files for implmentaion of pARMS
ISRCINC	                = -I../../src/include	

# path of the header files
IFFLAGS			=  -I../../include

# pARMS lib and its directory
#LIBDIR			=	$../../lib
#LIB			=	${LIBDIR}/libparms.a
#LIBFLAGS		=	-L${LIBDIR}
#SHFLAGS			=	
PARMS_LIBS		=	-L../../lib -lparms


allexe: dd-grid.ex fdd-grid.ex

dd-grid.ex: dd-grid.o skit.o functs.o fdmat.o aux.o
	${LINKER} ${CLINKFLAGS} -o $@ dd-grid.o skit.o functs.o fdmat.o aux.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} ${XLIBDIR} ${XLIB} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR}  ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} 

fdd-grid.ex: fdd-grid.o skit.o functs.o fdmat.o aux.o
	${F90LINKER} ${FLINKFLAGS} -o $@ fdd-grid.o skit.o functs.o fdmat.o aux.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} ${XLIBDIR} ${XLIB} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR}  ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} 

.c.o:
	${CC} ${ICFLAGS} ${ISRCINC} ${XIFLAGS} $(COPTFLAGS) \
	${CFLAGS} ${CFFLAGS} $< -c -o $@

.F.o:
	${F90} ${IFFLAGS} ${FFLAGS} $< -c -o $(@F)

.f.o:
	${FC} ${FFLAGS} $< -c -o $(@F)

clean:
	${RM} -rf *.o *.d core *~ *# *.ex paramfile.*




