SIESTA_ARCH=g95-nolibs
#
# Optimization options have to be investigated further
#
FC=g95
FC_ASIS=$(FC)
RANLIB=echo
#
FFLAGS= -O -Wall
FFLAGS_DEBUG= -g -O0 -Wall
LDFLAGS=
COMP_LIBS=linalg.a
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS=
SYS=
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#








