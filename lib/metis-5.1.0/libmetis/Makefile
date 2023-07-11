include ../../../Makefile.in
include ../Makefile.in


INCLUDES = -I. -I../GKlib -I../include

CFLAGS = $(COPTIONS) $(OPT) $(INCLUDES)

# OBJS= auxapi.o     coarsen.o     fm.o          kwayrefine.o  options.o  \
#      balance.o    compress.o    fortran.o     graph.o       minconn.o  \
#      parmetis.o   separator.o   bucketsort.o  contig.o      frename.o  \
#      initpart.o   mcutil.o      mincover.o    pmetis.o      sfm.o      \
#      timing.o     checkgraph.o  debug.o       gklib.o       kmetis.o   \
#      mesh.o       mmd.o         srefine.o     util.o        kwayfm.o   \
#      meshpart.o   ometis.o      refine.o      stat.o        wspace.o
OBJS=auxapi.o \
balance.o \
bucketsort.o \
checkgraph.o \
coarsen.o \
compress.o \
contig.o \
debug.o \
fm.o \
fortran.o \
frename.o \
gklib.o \
graph.o \
initpart.o \
kmetis.o \
kwayfm.o \
kwayrefine.o \
mcutil.o \
mesh.o \
meshpart.o \
minconn.o \
mincover.o \
mmd.o \
ometis.o \
options.o \
parmetis.o \
pmetis.o \
refine.o \
separator.o \
sfm.o \
srefine.o \
stat.o \
timing.o \
util.o \
wspace.o 


.c.o:
	$(CC) $(CFLAGS) -c $*.c

../lib/libmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

