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
#
SIESTA_ARCH=sgi64-mpich-cdf
FC=f90 -64
RANLIB=echo
#
FFLAGS=  -O3 -OPT:Olimit=0
NOOPT= 
FFLAGS_DEBUG= -g -O0
#
NETCDF_LIBS=-L/usr/local/netcdf-3.5/lib -lnetcdf
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/mpich-1.2/sgi64/include
MPI_LIBS= /usr/local/mpich-1.2/sgi64/lib/libmpich.a
BLACS_LIBS=-L/usr/local/lib/64 -lblacs
SCALAPACK_LIBS=-L/usr/local/lib/64 -lscalapack -lpblas -ltools -lredist
DEFS_MPI=-DMPI
#
LIBS= $(SCALAPACK_LIBS) $(BLACS_LIBS)  -lcomplib.sgimath  \
      $(MPI_LIBS) $(NETCDF_LIBS)
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $<
#




