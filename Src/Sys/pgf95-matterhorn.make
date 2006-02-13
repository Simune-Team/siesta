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
SIESTA_ARCH=pgf95-matterhorn
#
FC=pgf95
FC_ASIS=$(FC)
#
FFLAGS= -fast
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=         #  /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=    #  libnetcdf_f90.a
DEFS_CDF=            #  -DCDF
#
MPI_ROOT=/opt/64/mpich-gm-pgi6
MYRINET=/opt/64/gm-2.1.5/lib
MPI_LIBS= -L$(MPI_ROOT)/lib -lmpich -L$(MYRINET) -lgm
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/opt/64/mpich-gm-pgi6/include
DEFS_MPI=-DMPI
#

LIBS= -L/opt/64/acml-2.5.0/pgi64/lib  -lscalapack \
       -lblacs -lblacsF77init -lblacsCinit -lblacsF77init -lblacs \
       -llapack -lblas \
       $(MPI_LIBS)  $(NETCDF_LIBS)
SYS=cpu_time
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
# Compile atom.f and electrostatic.f without optimization.
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG) atom.f
#
electrostatic.o:
	$(FC) -c $(FFLAGS_DEBUG) electrostatic.f
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
