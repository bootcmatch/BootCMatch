# This is Make.inc file generated automatically by the PSBLAS configure script.
# It should be ready to use, included by Makefile.
# If it gives problems, consider editing it.

# These lines are quite portable.
.mod=.mod
.fh=.fh
.SUFFIXES: 
.SUFFIXES: .f90 .F90 .f .F .c .o


# The following ones are the variables used by the PSBLAS make scripts.

F90=gfortran 
FC=gfortran 
CC=gcc
F90COPT=-O3 -pg
#-pg 
FCOPT=-O3   -pg
#-pg
CCOPT=-O3   -pg
#-pg

LINKER=$(CC)
LIBS=-lsuperlu -L${SUPERLU_LIBDIR} -lhsl_mc64 -L${HSL_LIBDIR} -llapack -lf77blas -latlas  -L${HOME}/NUMERICAL/LIB/atlas/gnu491 -lgfortran
CINCLUDES=-DHAVE_SUPERLU -DSLU_VERSION_5 -I${SUPERLU_INCDIR} -DHAVE_HSL -I${HSL_INCDIR}

# These should be portable rules, arent they?
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<
.f.o:
	$(FC) $(FCOPT) $(FINCLUDES)  -c $<
.f90.o:
	$(F90) $(F90COPT) $(FINCLUDES) -c $<
.F.o:
	$(FC)  $(FCOPT)   $(FINCLUDES) $(FDEFINES) -c $<
.F90.o:
	$(F90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<

