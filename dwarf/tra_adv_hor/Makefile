# Horizontal tracer advection dwarf. FESOM's gfortran Release flags (without OpenMP: pure MPI).
FC      = mpif90
FFLAGS  ?= -march=native -O3 -finit-local-zero -finline-functions -fimplicit-none -ffree-line-length-none -cpp
EXTRA   ?=
OBJ = dwarf_data.o kernels_ref.o kernels_opt.o dwarf_main.o

dwarf.x: $(OBJ)
	$(FC) $(FFLAGS) $(EXTRA) -o $@ $(OBJ)

%.o: %.F90
	$(FC) $(FFLAGS) $(EXTRA) -c $<

kernels_ref.o kernels_opt.o: dwarf_data.o
dwarf_main.o: kernels_ref.o kernels_opt.o

vecreport:
	$(FC) $(FFLAGS) $(EXTRA) -fopt-info-vec-all=vec_ref.txt -c kernels_ref.F90 -o /dev/null
	$(FC) $(FFLAGS) $(EXTRA) -fopt-info-vec-all=vec_opt.txt -c kernels_opt.F90 -o /dev/null

clean:
	rm -f *.o *.mod dwarf.x
