.SUFFIXES:
.SUFFIXES: .F90 .F .f .o

CPP=cpp
FC=mpif90
LD=$(FC)
FCFLAGSFIXED=-g -c -O3 -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fconvert=swap -fopenmp $(NETCDF_INCLUDE) $(GRIB_API_INCLUDE)
FCFLAGSFREE=$(FCFLAGSFIXED)
CPPFLAGS=-traditional -P
LDFLAGS=-g -O3 -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fconvert=swap -fopenmp $(MAGPLUSLIB_SHARED) $(NETCDF_LIB) $(GRIB_API_LIB)
AR=ar
ARFLAGS=-rv

OBJ=ifs_modules.o ifs_interface.o ifs_notused.o

all: libfesom.a

.F90.o:
	$(CPP) $(CPPFLAGS) $< > $*.pp.f90
	$(FC) $(FCFLAGSFREE) $*.pp.f90 -o $*.o

.F.o:
	$(CPP) $(CPPFLAGS) $< > $*.pp.f
	$(FC) $(FCFLAGSFIXED) $*.pp.f -o $*.o

.f.o:
	$(CPP) $(CPPFLAGS) $< > $*.pp.f
	$(FC) $(FCFLAGSFIXED) $*.pp.f -o $*.o

libfesom.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)

clean:
	rm -f *.o *.mod *~ *.x *.pp.f90 *.pp.f *.lst *.in *.nc *.grb *.a
	rm -rf *.dSYM

ifs_interface.o: ifs_modules.o 
ifs_notused.o: ifs_modules.o


