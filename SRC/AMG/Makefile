#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
#
TOP=../..
include $(TOP)/make.inc

HEADERS =\
 bcm_amg.h

FILES =\
 bcm_amg.c \
 bcm_linmatch.c \
 bcm_weightmat.c \
 bcm_seqsuitor.c \
 bcm_adaptivecoarsening.c


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
