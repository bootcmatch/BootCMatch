#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
#
TOP=../..
include $(TOP)/make.inc

HEADERS =\
 bcm_linsolvers.h \
 bcm.h

FILES =\
 bcm_boot_solver.c \
 bcm_bootpcg.o

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
