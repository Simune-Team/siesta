# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996- .
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
# Makefile for Client-Server example
#
# JM Soler, J Nara, A Garcia, E Anglada
#
WORK_DIR=$(shell pwd)
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
MAIN_OBJDIR:=$(shell cd -P -- ../../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
TOP_LEVEL=$(shell pwd)../../../../../Src
VPATH:=$(TOP_LEVEL)
#
PSML_INCFLAGS=-I$(MAIN_OBJDIR)/psml/src
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(INCFLAGS)
#
default: protoNEB clean 
#
# Routines used by pipes version (VPATH knows where they are)
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
XC=libSiestaXC.a
#
include $(MAIN_OBJDIR)/Mk/FoX.mk
#
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a
PSML=$(MAIN_OBJDIR)/psml/src/libpsml.a
#
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	$(NCPS) $(PSML) $(PSOP) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(FOX_LIBS) $(LIBS)
#
libs_collected:
	(cd $(MAIN_OBJDIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WORK_DIR) )
	touch libs_collected
#
protoNEB: libs_collected protoNEB.o 
	$(FC) $(LDFLAGS) -o protoNEB \
	      protoNEB.o $(ALL_LIBS)
#
clean: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make pristine to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f protoNEB
	rm -f *.a *.o *.mod libs_collected

