.SUFFIXES:
.SUFFIXES: .F90 .F .f .o

CPP=cpp
FC=mpif90
LD=$(FC)
FCFLAGSFIXED=-g -c -O3 -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fconvert=swap -fopenmp $(NETCDF_INCLUDE) $(GRIB_API_INCLUDE)
FCFLAGSFREE=$(FCFLAGSFIXED)
CPPFLAGS=-traditional -P -Dkey_mpp_mpi
LDFLAGS=-g -O3 -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fconvert=swap -fopenmp $(MAGPLUSLIB_SHARED) $(NETCDF_LIB) $(GRIB_API_LIB)
AR=ar
ARFLAGS=-rv

OBJ=scripremap.o scripgrid.o parinter.o nctools.o ifs_modules.o ifs_interface.o ifs_notused.o

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

scripremap.o: nctools.o scrippar.o scripgrid.o
scripgrid.o: nctools.o scrippar.o
parinter.o: scripremap.o scrippar.o
interinfo.o: parinter.o
nemogcmcoup_mlflds_get.o: par_kind.o
nemogcmcoup_exflds_get.o: par_kind.o
nemogcmcoup_coup.o: par_kind.o
nemogcmcoup_coupinit.o: scripremap.o parinter.o interinfo.o
nemogcmcoup_wam_update.o: par_kind.o
nemogcmcoup_wam_get.o: par_kind.o
nemogcmcoup_wam_update_stress.o: par_kind.o
nemogcmcoup_mlinit.o: par_kind.o
nemogcmcoup_update_add.o: par_kind.o
nemogcmcoup_update.o: par_kind.o
nemogcmcoup_lim2_update.o: par_kind.o
nemogcmcoup_get.o: par_kind.o
nemogcmcoup_lim2_get.o: par_kind.o
