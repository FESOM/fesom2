include ../Makefile.in

LIB = ../lib/libmetis.a

INCLUDES = -I. -I../include

CFLAGS = $(COPTIONS) $(OPT) $(INCLUDES)

# OBJS= b64.o       csr.o    evaluate.o fs.o         gkregex.o htable.o  \
#      itemsets.o  memory.o pdb.o      random.o     seq.o     string.o  \
#      tokenizer.o blas.o   error.o    fkvkselect.o getopt.o  graph.o   \
#      io.o        mcore.o  omp.o      pqueue.o     rw.o      sort.o    \
#      timers.o    util.o
OBJS=b64.o \
blas.o \
csr.o \
error.o \
evaluate.o \
fkvkselect.o \
fs.o \
getopt.o \
gkregex.o \
graph.o \
htable.o \
io.o \
itemsets.o \
mcore.o \
memory.o \
omp.o \
pdb.o \
pqueue.o \
random.o \
rw.o \
seq.o \
sort.o \
string.o \
timers.o \
tokenizer.o \
util.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

../lib/libmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

