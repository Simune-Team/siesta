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

# Compiler and system module
# Create your own for a different compiler.
#
.SUFFIXES: .f90 .o

SYSMOD=

#  IBM
#FC=xlf90 -qsuffix=f=f90:cpp=F90
#SYSMOD=ibm_m_system.o

#  NAG
# FC=f95 -u
# SYSMOD=nag_m_system.o

#  Lahey
#FC=lf95
#SYSMOD=lahey_m_system.o

#FC=f95
#SYSMOD=alpha_m_system.o

#FC=pgf90
#SYSMOD=pgi_m_system.o

#  Intel
FC=ifc -w -mp -tpp5 -O3 -Vaxlib

default: driver

driver.o: fmixmd.o fsiesta.o sample.o
fmixmd.o: fsiesta.o sample.o
#

FFLAGS= -g -O0
FFLAGS_DEBUG= -g
LDFLAGS=

OBJS= $(SYSMOD) fmixmd.o fsiesta.o sample.o 

# Macros
.f90.o:
	$(FC) -c $(FFLAGS) $<

driver: $(OBJS) driver.o
	$(FC) -o driver driver.o $(OBJS)

siesta:
	(cd ../../../../Src ; $(MAKE) siesta)

clean: 	
	rm -f $(OBJS) *.mod
	rm -f driver *.o
