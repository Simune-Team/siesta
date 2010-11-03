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
SIESTA_DIR=$(WORK_DIR)/../../../..
SRC_DIR=$(SIESTA_DIR)/Src
OBJ_DIR=$(SIESTA_DIR)/Obj-gfortran-parallel
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH=$(SRC_DIR)
#
#  simple_pipes_serial   : Simple serial test with pipes version
#  simple_pipes_parallel : Simple parallel test with pipes version
#  simple_mpi_serial     : Simple serial test with MPI version
#  simple_mpi_parallel   : Simple parallel test with MPI version
#
default: simple_pipes_serial simple_pipes_parallel \
	simple_mpi_serial simple_mpi_parallel clean
#
ARCH_MAKE=$(OBJ_DIR)/arch.make
include $(ARCH_MAKE)
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
# Routines used by pipes version (VPATH knows where they are)
#
FSIESTA= fsiesta_pipes.o pxf.o
#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
XC=libSiestaXC.a
FoX_LIBS=`$(OBJ_DIR)/FoX/FoX-config --libs --wcml`
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS) $(LIBS)
#
libs_collected:
	(cd $(OBJ_DIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WORK_DIR) )
	touch libs_collected
#
# Pipes version
#
simple_pipes_serial: $(FSIESTA) simple_serial.o
	$(FC) $(LDFLAGS) -o simple_pipes_serial $(FSIESTA) simple_serial.o 
#
simple_pipes_parallel: $(FSIESTA) simple_pipes_parallel.o
	$(FC) $(LDFLAGS) -o simple_pipes_parallel \
	      $(FSIESTA) simple_pipes_parallel.o
#
# MPI version
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
	rm -f simple_mpi_serial simple_mpi_parallel
	rm -f *.a *.o *.mod libs_collected

