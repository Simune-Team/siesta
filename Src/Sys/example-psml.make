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
SIESTA_ARCH=gfortran-nolibs-netcdf
#
#
FC=gfortran
#
FC_ASIS=$(FC)
#
FFLAGS=-g -fbacktrace -O2
#FFLAGS= -g -O0 -Wall -fcheck=all -fbacktrace
FFLAGS_DEBUG= -g -Wall -Wextra
LDFLAGS=
RANLIB=echo
LIBS=  
SYS=nag
#
FOX_ROOT=$(HOME)/lib/FoX/gfortran
#
# --- Edit the location of your netcdf files
#
NETCDF_ROOT=/usr/local
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
#
PSML_ROOT=$(HOME)/S/GIT/psml
PSML_INCFLAGS=-I$(PSML_ROOT)/src
#
FPPFLAGS=-DGFORTRAN -DCDF -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 
COMP_LIBS=
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
PSML_LIBS= -L$(PSML_ROOT)/src -lpsml
LIBS=$(NETCDF_LIBS) $(PSML_LIBS) -framework veclib
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#

