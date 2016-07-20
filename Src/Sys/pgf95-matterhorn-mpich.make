# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
SIESTA_ARCH=pgf95-matterhorn-mpich
#
FC=pgf95
FC_ASIS=$(FC)
#
FFLAGS= -fastsse
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
FPPFLAGS_CDF=            #  -DCDF
#
MPI_ROOT_MPICH=/opt/64/mpich-pgi6
## Be sure to use the right mpirun later on...
MPIRUN=/opt/64/mpich-pgi6/bin/mpirun
#
# To use myrinet libraries instead of default mpich ones use the symbols below
#MPI_ROOT_GM=/opt/64/mpich-gm-pgi6
#MYRINET=/opt/64/gm-2.1.5/lib
#MYRINET_LIBS=-L$(MYRINET) -lgm
#
MPI_ROOT=$(MPI_ROOT_MPICH)
MPI_LIBS= -L$(MPI_ROOT)/lib -lmpich 
##MPI_LIBS= -L$(MPI_ROOT)/lib -lmpich $(MYRINET_LIBS)
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=$(MPI_ROOT)/include
FPPFLAGS_MPI=-DMPI
#

LIBS= -L/opt/64/acml-2.5.0/pgi64/lib  -lscalapack \
       -lblacs -lblacsF77init -lblacsCinit -lblacsF77init -lblacs \
       -llapack -lblas \
       $(MPI_LIBS)  $(NETCDF_LIBS)
SYS=cpu_time
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
# Compile atom.f and electrostatic.f without optimization.
# Make sure to use an explicit dependency and $<, so that these
# lines work with VPATH.
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
