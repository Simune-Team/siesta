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
SIESTA_ARCH=g95-macosx-netcdf
#
#
FC=g95
#
FC_ASIS=$(FC)
#
FFLAGS=-O0 -g 
FFLAGS_DEBUG= -g -Wall
LDFLAGS=
RANLIB=echo
METIS_LIBS=/Users/ag/lib/metis/lib/libmetis.a
SYS=nag
#
# --- Edit the location of your netcdf files
#
NETCDF_ROOT=$(HOME)/lib/netcdf-3.6.2/g95
INCFLAGS=-I$(NETCDF_ROOT)/include
#
DEFS=-DCDF -DFC_HAVE_FLUSH -DFC_HAVE_ABORT  -DDEBUG   # Note this !!
COMP_LIBS=linalg.a
#
NETCDF_LIBS= -L$(NETCDF_ROOT)/lib -lnetcdf
LIBS=$(METIS_LIBS) $(NETCDF_LIBS) -framework veclib
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#

