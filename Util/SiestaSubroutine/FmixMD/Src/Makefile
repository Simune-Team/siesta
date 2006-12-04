# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
# Makefile for Client-Server example
#
# JM Soler, A Garcia
#
.SUFFIXES:
.SUFFIXES: .f90 .o 
#
#       driver: Mixed molecular dynamics
#       simple: Simple test
#       para  : Test driver for parallel calculation
#
default: what simple driver para
#
# You should have compiled Siesta in the main Src directory before
# attempting to compile these examples.
#
include ../../../../Src/arch.make
#
# Uncomment the following line for debugging support
#
FFLAGS=$(FFLAGS_DEBUG)
#
FSIESTA_OBJ=../../../../Src/fsiesta.o      # Client helper routines
FSIESTA_MOD=../../../../Src/fsiesta.mod    # Client helper module
PXF_OBJ=../../../../Src/pxf.o              # Portable flush

#
what:
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_ARCH)"
	@echo
	cp -f $(FSIESTA_OBJ) $(FSIESTA_MOD) $(PXF_OBJ) .
#
#
driver.o: fmixmd.o sample.o
fmixmd.o: sample.o

#
OBJS=sample.o fmixmd.o driver.o fsiesta.o pxf.o
driver : $(OBJS)
	$(FC) -o driver \
              $(LDFLAGS) $(OBJS) 
#
simple: simple.o fsiesta.o pxf.o
	$(FC) $(LDFLAGS) -o simple simple.o fsiesta.o pxf.o
#
para: para.o fsiesta.o pxf.o
	$(FC) $(LDFLAGS) -o para para.o fsiesta.o pxf.o
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f driver simple para *.o *.mod



