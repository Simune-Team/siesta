# ---
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone unfold
#

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

default: unfold

WDIR=$(shell pwd)
VPATH=../../../Src
OBJDIR=../../../Obj
HSXDIR=../../HSX
ARCH_MAKE=$(OBJDIR)/arch.make
include $(ARCH_MAKE)

SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
WXML=libwxml.a
XC=libSiestaXC.a
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	  $(MPI_INTERFACE) $(COMP_LIBS) $(LIBS)

libs_collected:
	(cd $(OBJDIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WDIR) )
	cp -f $(HSXDIR)/hsx_m.f90 $(WDIR)
	touch libs_collected

unfold.o: unfold.F90 hsx_m.f90 libs_collected
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F) hsx_m.f90 unfold.F90

unfold: libs_collected unfold.o
	$(FC) $(LDFLAGS) -o unfold hsx_m.o unfold.o $(ALL_LIBS)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f unfold *.o  *.a *.mod libs_collected

