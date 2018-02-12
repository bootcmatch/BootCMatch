#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
#
TOP=../..
include $(TOP)/make.inc

HEADERS =\
 bcm_matvec.h 

FILES =\
 bcm_csr_matop.c\
 bcm_csr_matrix.c\
 bcm_csr_matvec.c\
 bcm_csr_relax.c\
 bcm_vector.c 

INCDIR=$(TOP)/include
CINCLUDES += -I$(INCDIR)
OBJS = $(FILES:.c=.o)

##################################################################
# Targets
##################################################################

all: lib copy

copy:
	/bin/cp -p $(HEADERS) $(INCDIR)

clean:
	rm -f *.o 

##################################################################
# Rules
##################################################################

lib: $(OBJS)
	ar cur $(TOP)/$(LIBNAME) $(OBJS)
	ranlib $(TOP)/$(LIBNAME)

$(OBJS) $(FILES): $(HEADERS)
