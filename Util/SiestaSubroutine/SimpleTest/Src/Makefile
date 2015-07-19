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
# JM Soler, A Garcia, E Anglada
#
WORK_DIR=$(shell pwd)
OBJDIR=Obj
#
.SUFFIXES:
.SUFFIXES: .c .f .F .o .a .f90 .F90
#
MAIN_OBJDIR:=$(shell cd -P -- ../../../../\$(OBJDIR) && pwd -P)
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
#
TOP_LEVEL=$(shell pwd)/../../../../Src
VPATH:=$(TOP_LEVEL)
#
PSML_INCFLAGS=-I$(MAIN_OBJDIR)/psml/src
NCPS_INCFLAGS=-I$(MAIN_OBJDIR)/ncps/src
PSOP_INCFLAGS=-I$(MAIN_OBJDIR)/psoplib/src
INCFLAGS:=$(NETCDF_INCFLAGS) \
          $(NCPS_INCFLAGS) $(PSOP_INCFLAGS) $(PSML_INCFLAGS) \
          $(INCFLAGS)
#
#  simple_pipes_serial     : Simple serial test with pipes version
#  simple_pipes_parallel   : Simple parallel test with pipes version
#  simple_sockets_serial   : Simple serial test with sockets version
#  simple_sockets_parallel : Simple parallel test with sockets version
#  simple_mpi_serial       : Simple serial test with MPI version
#  simple_mpi_parallel     : Simple parallel test with MPI version
#
# default: simple_pipes_serial simple_sockets_serial clean

default: simple_pipes_serial simple_pipes_parallel \
	simple_sockets_serial simple_sockets_parallel \
	simple_mpi_serial simple_mpi_parallel clean 
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
# Routines used by pipes and sockets versions (VPATH knows where they are)
#
FSIESTA_PIPES= fsiesta_pipes.o pxf.o
FSIESTA_SOCKETS= fsockets.o sockets.o fsiesta_sockets.o 
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
XC=libSiestaXC.a
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
# Pipes version
#
simple_pipes_serial: $(FSIESTA_PIPES) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_pipes_serial $(FSIESTA_PIPES) simple_serial.o 
#
simple_pipes_parallel: $(FSIESTA_PIPES) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_pipes_parallel \
	      $(FSIESTA_PIPES) simple_parallel.o
#
# Sockets version
#
simple_sockets_serial: $(FSIESTA_SOCKETS) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_sockets_serial \
	      $(FSIESTA_SOCKETS) simple_serial.o 
#
simple_sockets_parallel: $(FSIESTA_SOCKETS) simple_parallel.o
	$(FC) $(LDFLAGS) -o simple_sockets_parallel \
	      $(FSIESTA_SOCKETS) simple_parallel.o
#
# MPI version
#
phonons: libs_collected phonons.o 
	$(FC) $(LDFLAGS) -o phonons phonons.o $(ALL_LIBS)
#
simple_mpi_serial: libs_collected simple_serial.o 
	$(FC) $(LDFLAGS) -o simple_mpi_serial simple_serial.o $(ALL_LIBS)
#
simple_mpi_parallel: libs_collected simple_mpi_parallel.o 
	$(FC) $(LDFLAGS) -o simple_mpi_parallel \
	      simple_mpi_parallel.o $(ALL_LIBS)
#
clean: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make pristine to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f simple_pipes_serial simple_pipes_parallel 
	rm -f simple_sockets_serial simple_sockets_parallel 
	rm -f simple_mpi_serial simple_mpi_parallel
	rm -f phonons
	rm -f *.a *.o *.mod libs_collected

