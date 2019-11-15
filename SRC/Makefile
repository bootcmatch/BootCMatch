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


libs:  mutilities  mmatvec mamg mbootamg   mlinsolvers


###################################################################
# Rules
# ##################################################################
#

mbootamg:
	( cd BOOTAMG; $(MAKE) LIBNAME=$(LIBNAME))

mutilities:
	( cd utilities; $(MAKE)  LIBNAME=$(LIBNAME))

mmatvec:
	( cd MATVEC; $(MAKE)  LIBNAME=$(LIBNAME))

mamg:
	( cd AMG; $(MAKE)  LIBNAME=$(LIBNAME))

mlinsolvers:
	( cd LINSOLVERS; $(MAKE)  LIBNAME=$(LIBNAME))

clean: 
	( cd utilities; $(MAKE) clean) 
	( cd MATVEC; $(MAKE) clean) 
	( cd BOOTAMG;  $(MAKE) clean) 
	( cd AMG; $(MAKE) clean)
	( cd LINSOLVERS; $(MAKE) clean)
