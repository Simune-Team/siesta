SIESTA_ARCH=lahey-mpich
#
FC=lf95
FC_ASIS=$(FC)
FFLAGS= -O  --ntrace --tpp
LDFLAGS=--staticlink
FFLAGS_DEBUG= -g -O0 --chk
#
NETCDF_LIBS=-L/usr/local/netcdf-3.5/lib/lahey -lnetcdf
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/mpich-1.2/lahey/include
DEFS_MPI=-DMPI
#
LIBS= -L/usr/local/lib/lahey \
      -lscalapack -lpblas \
      -ltools -lredist \
      -lblacs.mpich  \
      -llapack -lblas \
      /usr/local/mpich-1.2/lahey/lib/libmpich.a $(NETCDF_LIBS)
SYS=bsd
DEFS= $(DEFS_MPI) $(DEFS_CDF)
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
