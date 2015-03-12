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
SIESTA_ARCH=gfortran-macosx-netcdf-psml
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
# This is an external version of FoX. If you
# uncomment this, the internal version will be
# used
#
FOX_ROOT=$(HOME)/lib/FoX/gfortran
#
# --- Edit the location of your netcdf files
#
NETCDF_ROOT=/usr/local
NETCDF_INCFLAGS=-I$(NETCDF_ROOT)/include
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdff
#
# --- Edit the location of your psml files
#
PSML_ROOT=$(HOME)/S/GIT/psml
PSML_INCFLAGS=-I$(PSML_ROOT)/src
PSML_LIBS= -L$(PSML_ROOT)/src -lpsml
#
FPPFLAGS=-DGFORTRAN -DCDF -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 
COMP_LIBS=
#
# These are the external libs only
#
LIBS=$(NETCDF_LIBS) $(PSML_LIBS)  -framework veclib
#
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

