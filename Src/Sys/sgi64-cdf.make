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
SIESTA_ARCH=sgi64-cdf
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
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS=  -lcomplib.sgimath $(NETCDF_LIBS)
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




