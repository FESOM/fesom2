FC = gfortran
CXX = g++
#todo: should be conditional since in the cluster this is already set and will overwrite the flag
HECUBA_ROOT = /usr/local/hecuba
FCFLAGS = -Wall -Wextra -cpp
#CCFLAGS = -Wall -Wextra -lstdc++
CCFLAGS = -Wall -Wextra -lstdc++ -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib
LDFLAGS = -lstdc++ -I${HECUBA_ROOT}/include -L${HECUBA_ROOT}/lib -lhfetch -Wl,-rpath,${HECUBA_ROOT}/lib

all: test.x
test.o : hecuba_exporter_mod.o

%.x : %.o hecuba_exporter_mod.o Hecuba_exporter_wrapper.o Hecuba_exporter.o
	${FC} $^ -o $@ ${LDFLAGS}

%.o : %.f90
	${FC} ${FCFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CCFLAGS} -c $^ -o $@

.PHONY : clean

clean :
	${RM} -rf *.o *.mod test.x
