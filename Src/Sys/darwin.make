# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=darwin
#
FC=f90
FC_ASIS=$(FC)
#
FFLAGS=  -O
FFLAGS_DEBUG= -g
RANLIB=echo "... do we need to ranlib this? :"
LDFLAGS=
COMP_LIBS=linalg.a
#
NETCDF_LIBS=
NETCDF_INTERFACE=
FPPFLAGS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
LIBS= -lblas -lU77 $(NETCDF_LIBS)
SYS=bsd
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
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
