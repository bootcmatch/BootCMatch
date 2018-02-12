#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
#
TOP=../..
include $(TOP)/make.inc

HEADERS =\
 bcm_bootamg.h

FILES =\
 bcm_bootamg.c \
 bcm_bootstrap.c \
 bcm_boot_prec.c\
 bcm_inneritkcycle.c\
 bcm_kc_apply.c  \
 bcm_gcycle.c


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
