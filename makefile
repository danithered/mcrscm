
IDIR =./src/include
ODIR=./src/obj
SRCDIR=./src

CC=g++ -std=c++17
C=gcc

CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` `pkg-config --cflags RNAlib2` -O3 -lboost_system -lboost_serialization  # for stuff

LIBS=-lm `pkg-config --libs gsl` `pkg-config --libs RNAlib2` -fopenmp # for system RNAlib

model_mcrs=-D "_MCRS_"
model_scm=-D "_SCM_"

# receipts

_DEPS_COMMON = randomgen.h dv_tools.h rnarep.h annot.h parameters.h

_DEPS_MCRS = ca.h randomgen.h dv_tools.h mcrs_parameters.h rnarep.h annot.h parameters.h
_OBJ_MCRS = mcrs.o ca.o dv_tools.o mcrs_parameters.o rnarep.o annot.o randomgen.o

_DEPS_SCM = broken.hpp randomgen.h dv_tools.h parameters.h scm_parameters.h rnarep.h annot.h cm.h rnarep_serialise.h cm_serialise.h  
_OBJ_SCM = scm.o dv_tools.o scm_parameters.o rnarep.o annot.o randomgen.o cm.o 

_OBJ_randseq = randseq.o randomgen.o rnarep.o annot.o mcrs_parameters.o dv_tools.o

_OBJ_rev = reverse.o ca.o annot.o rnarep.o mcrs_parameters.o randomgen.o dv_tools.o 

_OBJ_strgen = strgen.o randomgen.o 

# assemble deps

DEPS_MCRS = $(patsubst %,$(IDIR)/%,$(_DEPS_MCRS))
OBJ_MCRS = $(patsubst %,$(ODIR)/%,$(_OBJ_MCRS))
DEPS_SCM = $(patsubst %,$(IDIR)/%,$(_DEPS_SCM))
OBJ_SCM = $(patsubst %,$(ODIR)/%,$(_OBJ_SCM))
OBJ_randseq = $(patsubst %,$(ODIR)/%,$(_OBJ_randseq))
OBJ_rev = $(patsubst %,$(ODIR)/%,$(_OBJ_rev))
OBJ_strgen = $(patsubst %,$(ODIR)/%,$(_OBJ_strgen))


# def values
DEPS=$(DEPS_MCRS)
model=$(model_mcrs)

# compile o files

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	@mkdir -p ${ODIR}
	$(CC) -c -o $@ $< $(CFLAGS) $(model)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p ${ODIR}
	$(C) -c -o $@ $< $(CFLAGS) $(model)

# executables

mcrs: $(OBJ_MCRS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: scm
scm: model=$(model_scm)
scm: DEPS=$(DEPS_SCM)
scm: $(OBJ_SCM)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: randseq
randseq: model=$(model_mcrs)
randseq: DEPS=$(DEPS_MCRS)
randseq: $(OBJ_randseq)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: rev
rev: $(OBJ_rev)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: gen
gen: $(OBJ_strgen)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

.PHONY: all
.DEFAULT_GOAL := all
all: mcrs randseq scm

