# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=cscs-cray
#
# For Cray XT-3 at CSCS (Palu)
# It uses pgf90 V6 (f95 !) as a backend
#
FC=ftn -target=catamount
FPP=linux-pgf90 -F
FC_ASIS=$(FC)
#
FFLAGS= -fast
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=            #  -DCDF
#
MPI_LIBS=
MPI_INTERFACE=
MPI_INCLUDE=
FPPFLAGS_MPI=
#
LIBS= 
SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
#
# Important 
# Compile atom.f and electrostatic.f without optimization.
# Make sure that the dependency is explicit, so that these
# lines work with VPATH
#
atom.o: atom.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
electrostatic.o: electrostatic.f
	$(FC) -c $(FFLAGS_DEBUG) $<
#
.F.o:
	$(FPP) $(FPPFLAGS) $<  ; mv $*.f aux_$*.f
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) -o $*.o aux_$*.f
	rm -f aux_$*.f
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FPP) $(FPPFLAGS) $<  ; mv $*.f aux_$*.f90
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) -o $*.o aux_$*.f90
	rm -f aux_$*.f90
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
