SHELL	= /bin/sh

include makefile.in

allexe: test.ex ftest.ex dd-petsc.ex

test.ex: test.o  parms.o chkopts
	${LINKER} ${CLINKFLAGS} -o $@ test.o  parms.o  \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} ${XLIBDIR} ${XLIB} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR}  ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} ${PETSC_KSP_LIB} 
	${RM} test.o

ftest.ex: ftest.o parms.o fparms.o chkopts
	${F90LINKER} ${FLINKFLAGS} -o $@ ftest.o parms.o fparms.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} ${XLIBDIR} ${XLIB} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR}  ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} ${PETSC_KSP_LIB}
	${RM} ftest.o 

dd-petsc.ex: dd-petsc.o parms.o readmat.o skit.o chkopts
	${LINKER} ${CLINKFLAGS} -o $@ dd-petsc.o parms.o readmat.o skit.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} ${XLIBDIR} ${XLIB} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR}  ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} ${PETSC_KSP_LIB}
	${RM} dd-petsc.o

cleanall::
	  @-${RM} -rf *.o core *~ *# *.ex paramfile.*


