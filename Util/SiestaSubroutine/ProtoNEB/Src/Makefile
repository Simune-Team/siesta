# Siesta/Util/SiestaSubroutine/SimpleTest/Makefile
# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
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

all: default
default: protoNEB

MAIN_OBJDIR:=$(shell cd -P -- ../../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
TOP_LEVEL=$(shell pwd)/../../../../Src
VPATH:=$(TOP_LEVEL)
#
#
# Routines used by pipes version (VPATH knows where they are)
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
#
MS=MatrixSwitch.a
XMLPARSER=$(XMLF90_LIBS)
XC=$(GRIDXC_LIBS)
PSML=$(PSML_LIBS)
#
NCPS=$(MAIN_OBJDIR)/ncps/src/libncps.a
PSOP=$(MAIN_OBJDIR)/psoplib/src/libpsop.a
#
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
#
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(XMLF90_INCFLAGS) $(GRIDXC_INCFLAGS) $(INCFLAGS)
#
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(XMLF90_LIBS)  $(XC) \
	$(NCPS) $(PSML) $(PSOP) $(MS)\
	$(MPI_INTERFACE) $(COMP_LIBS) $(LIBS)
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
intermediate: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make clean to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: clean
clean:
	@echo "==> Cleaning all intermediate and executable files"
	rm -f protoNEB
	rm -f *.a *.o *.mod libs_collected

