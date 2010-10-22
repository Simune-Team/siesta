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
OBJ_DIR=$(SIESTA_DIR)/Obj_parallel
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH=$(SRC_DIR)
#
#  pipes_serial   : Simple serial test with pipes version
#  pipes_parallel : Simple parallel test with pipes version
#  mpi_serial     : Simple serial test with MPI version
#  mpi_parallel   : Simple parallel test with MPI version
#
default: pipes_serial pipes_parallel mpi_serial mpi_parallel
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
pipes_serial: $(FSIESTA) serial.o
	$(FC) $(LDFLAGS) -o pipes_serial $(FSIESTA) serial.o 
#
pipes_parallel: $(FSIESTA) pipes_parallel.o
	$(FC) $(LDFLAGS) -o pipes_parallel $(FSIESTA) pipes_parallel.o
#
# MPI version
#
mpi_serial: libs_collected serial.o 
	$(FC) $(LDFLAGS) -o mpi_serial serial.o $(ALL_LIBS)
#
mpi_parallel: libs_collected mpi_parallel.o 
	$(FC) $(LDFLAGS) -o mpi_parallel mpi_parallel.o $(ALL_LIBS)
#
clean: 
	@echo "==> Cleaning intermediate files ONLY"
	@echo "==> Use make pristine to remove executable files"
	rm -f *.a *.o *.mod libs_collected
#
pristine: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f pipes_serial pipes_parallel mpi_serial mpi_parallel *.a *.o *.mod libs_collected

