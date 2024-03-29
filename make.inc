# This is Make.inc file.
# Eventually it will be generated with a configure script.

# 
.SUFFIXES: .f90 .F90 .f .F .c .o


#
# Optional packages
#

# SuperLU: sparse factorization
#  Versions 4 and 5 have a slightly different interface, need
#  to define which one. 
#
#SUPERLU_DIR=/opt/superlu/5.0/gnu/5.4.0
#SUPERLU_INCDIR=$(SUPERLU_DIR)/include
#SUPERLU_LIBDIR=$(SUPERLU_DIR)/lib
#SUPERLU_LIBS=-lsuperlu -L$(SUPERLU_LIBDIR)
#SUPERLU_LIBS=-lsuperlu
#SUPERLU_INCDIR=/usr/include/SuperLU
#SUPERLU_FLAGS=-DHAVE_SUPERLU -DSLU_VERSION_5 -I$(SUPERLU_INCDIR)

# HSL library (for MC64)
#HSL_DIR=/opt/hsl
#HSL_DIR=/opt/hsl/2.3.1/sys
#HSL_INCDIR=$(HSL_DIR)/include
#HSL_LIBDIR=$(HSL_DIR)/lib
#HSL_LIBS=-lhsl_mc64 -L$(HSL_LIBDIR)
#HSL_FLAGS= -DHAVE_HSL -I$(HSL_INCDIR)

# SPRAL package for auction algorithm
#SPRAL_DIR=/opt/spral
#SPRAL_DIR=/opt/spral/2015.04.20/sys
#SPRAL_INCDIR=$(SPRAL_DIR)/include
#SPRAL_LIBDIR=$(SPRAL_DIR)/lib
#SPRAL_LIBS=-lspral -L$(SPRAL_LIBDIR)
#SPRAL_FLAGS=-DHAVE_SPRAL -I$(SPRAL_INCDIR)


#
# LAPACK and BLAS are mandatory if you use SuperLU
#
#LAPACK_LIBS=-llapack
#BLAS_LIBS=-latlas -L/usr/lib64/atlas

#
# Compilers & such
#
FC=gfortran 
CC=gcc
F90COPT=-O3 
#-pg 
FCOPT=-O3  
#-pg
CCOPT=-O3  
#-pg

LINKER=$(CC) -pg
LIBS=$(SUPERLU_LIBS) $(HSL_LIBS) $(SPRAL_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

CINCLUDES=$(SUPERLU_FLAGS) $(HSL_FLAGS) $(SPRAL_FLAGS)

# These should be portable rules, arent they?
.c.o:
	$(CC) $(CCOPT) $(CFLAGS) $(CINCLUDES) -c $< -o $@
.f.o:
	$(FC) $(FCOPT) $(FINCLUDES)  -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FINCLUDES) -c $< -o $@
.F.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@

