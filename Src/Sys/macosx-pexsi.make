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
SIESTA_ARCH=gfortran-macosx64-openmpi
# The only thing you should change is the location of the libraries
# on your computer
#
FC=mpif90
#
FC_ASIS=$(FC)
#
FFLAGS= -g -O0 -fbacktrace #-fcheck=all
FFLAGS_CHECKS= -O0 -g -fcheck=all
FFLAGS_DEBUG= -g -O0
RANLIB=echo
COMP_LIBS=
#
#FOX_ROOT=$(HOME)/lib/FoX/4.1.2/gfortran-4.8.3
#
NETCDF_INCFLAGS=-I$(NETCDF_INCLUDE)
NETCDF_LIBS= $(NETCDF_FORTRAN_LIBS)
FPPFLAGS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=.      # Note . for no-op
FPPFLAGS_MPI=-DMPI

PEXSI_LIB=$(HOME)/lib/PEXSI/0.7.2/openmpi-1.8.1-gfortran-4.8.3/lib/libpexsi_osx_v0.7.2.a
SUPERLU_LIB=$(HOME)/lib/SuperLU_DIST_3.3/lib/libsuperlu_dist_3.3.a
PARMETIS_LIB=-L$(HOME)/lib/parmetis-4.0.2/lib -lparmetis
METIS_LIB=-L$(HOME)/lib/metis/lib -lmetis
MPICXX_LIB=/opt/openmpi-1.8.1-gfortran-4.8.3/lib/libmpi_cxx.dylib  -lstdc++ 

LIBS=-L/opt/scalapack/openmpi-1.6.1-gfortran/lib \
        -lscalapack -ltmg -lreflapack -lrefblas \
      $(NETCDF_LIBS) \
      $(PEXSI_LIB) $(SUPERLU_LIB) \
      $(PARMETIS_LIB) $(METIS_LIB) $(MPICXX_LIB)

SYS=nag
FPPFLAGS= $(FPPFLAGS_CDF) $(FPPFLAGS_MPI) -DMPI_TIMING
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
