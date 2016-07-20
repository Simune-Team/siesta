# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=atto-pgf95-openmpi
#
# PGI compiler, with OpenMPI and home-compiled BLACS and SCALAPACK
#
# The only thing you should change is the location of the PGI libraries
# on your computer
#
FC=/share/apps/openmpi-1.4.3-pgi/bin/mpif90
FC_ASIS=$(FC)
#
FFLAGS= -fast             # Other options might be useful
FFLAGS_CHECKS= -O0 -g -Mchkptr -Mbounds -traceback
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI
#
METIS_LIB=/share/apps/metis-4.0/libmetis.a

LAPACK=-L$(PGI)/osx86-64/10.8/lib -llapack -lblas
LIBS=-L/share/apps/local_libs/pgi-10.8 -lscalapack \
     -L/share/apps/local_libs/pgi-10.8/openmpi-1.4.2 -lblacsF77 -lblacsC -lblacs -lblacsF77 -lblacsC \
      $(LAPACK) \
      $(NETCDF_LIBS) $(METIS_LIB)
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
