#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
include ../make.inc
#
TOP=..
INCDIR=$(TOP)/include
LIBDIR=$(TOP)/lib
CINCLUDES +=-I$(INCDIR)

RUNDIR=./RUNS

all: test

test: testkbootsolvers testsolving testbootstrap

clean:
	rm -f *.o testkbootsolvers testbootstrap testsolving

##################################################################
# Rules
##################################################################
UTILS=ioutil.o

testkbootsolvers:  testkbootsolvers.o $(UTILS)
	$(LINKER) $(LINKOPT) -o $@ $< $(UTILS) -lBCM -L$(LIBDIR) $(LIBS) -lm $(LINKOPT)
	mv $@ $(RUNDIR)

testsolving:  testsolving.o $(UTILS)
	$(LINKER) $(LINKOPT) -o $@ $< $(UTILS) -lBCM -L$(LIBDIR) $(LIBS) -lm $(LINKOPT)
	mv $@ $(RUNDIR)

testbootstrap:  testbootstrap.o $(UTILS)
	$(LINKER) $(LINKOPT) -o $@ $<  $(UTILS) -lBCM -L$(LIBDIR) $(LIBS) -lm $(LINKOPT)
	mv $@ $(RUNDIR)
