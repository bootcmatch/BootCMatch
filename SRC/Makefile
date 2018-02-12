#BHEADER**********************************************************************
# # This file is part of BootCMatch.
# #
# # $release: 1.0 $
# #EHEADER**********************************************************************
TOP=..
include $(TOP)/make.inc
##################################################################
# Targets
# ##################################################################
LIBNAME=libBCM.a

all:  libs


libs:  util  matvec amg bootamg   linsolvers


###################################################################
# Rules
# ##################################################################
#

bootamg:
	( cd BOOTAMG; $(MAKE) LIBNAME=$(LIBNAME))

util:
	( cd utilities; $(MAKE)  LIBNAME=$(LIBNAME))

matvec:
	( cd MATVEC; $(MAKE)  LIBNAME=$(LIBNAME))

amg:
	( cd AMG; $(MAKE)  LIBNAME=$(LIBNAME))

linsolvers:
	( cd LINSOLVERS; $(MAKE)  LIBNAME=$(LIBNAME))

clean: 
	( cd utilities; $(MAKE) clean) 
	( cd MATVEC; $(MAKE) clean) 
	( cd BOOTAMG;  $(MAKE) clean) 
	( cd AMG; $(MAKE) clean)
	( cd LINSOLVERS; $(MAKE) clean)
