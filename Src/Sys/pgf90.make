SIESTA_ARCH=pgf90
#
FC=pgf90
FC_ASIS=$(FC)
#
FFLAGS=  -fast
FFLAGS_DEBUG= -g -O0
LDFLAGS=
COMP_LIBS=
#
NETCDF_LIBS=/usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= -L/usr/local/lib -llapack -lblas -lg2c $(NETCDF_LIBS)
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
