SIESTA_ARCH=pgf90-mpich
#
FC=ff90
FC_ASIS=$(FC)
#
FFLAGS= -fast
FFLAGS_DEBUG= -g -O0
#
NETCDF_LIBS=/usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
DEFS_MPI=-DMPI
#
LIBS= -L/usr/local/lib/pgi \
      -lscalapack -l1upblas -l1utools -l.pgi.aux -lredist \
      -lfblacs  -llapack -lblas \
       -l1umpich -lpgiarg $(NETCDF_LIBS)
SYS=bsd
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
