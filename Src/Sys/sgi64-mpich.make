#
SIESTA_ARCH=sgi64-mpich
FC=f90 -64
#
FFLAGS=  -O3 -OPT:Olimit=0
NOOPT= 
FFLAGS_DEBUG= -g -O0
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
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
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<
#




