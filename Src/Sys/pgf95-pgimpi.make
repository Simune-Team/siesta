# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=pgf95-pgimpi
#
# It uses PGI's own MPI subsystem
#
# The only thing you should change is the location of the PGI libraries
# on your computer
#
FC=pgf95 -Mmpi
FC_SERIAL=pgf95
#
FFLAGS= -fast             # Other options might be useful
#FFLAGS= -O0 -g -Mchkptr -Mbounds -traceback
FFLAGS_CHECKS= -O0 -g -Mchkptr -Mbounds -traceback
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI
#
# LAPACK is not needed if -Mscalapack is used ...
### LAPACK=-L$(PGI)/linux86-64/2010/lib -llapack -lblas
#
METIS_LIB=/share/apps/metis-4.0/libmetis.a
LIBS= -Mscalapack $(LAPACK) $(NETCDF_LIBS) $(METIS_LIB)
#
SYS=cpu_time
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
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
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)  $(FPPFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS)   $<
#
