SIESTA_ARCH=intel9-mkl8
#
# Intel fortran compiler 9 for linux with mkl 8 optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
# Note: The -mp1 option is necessary to recover IEEE floating point precision.
#
FC=ifort
#
FFLAGS= -w -xP -O3 -mp1
FFLAGS_DEBUG= -g 
LDFLAGS=
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
GUIDE=/opt/intel/mkl/8.0.1/lib/32/libguide.a
LAPACK=/opt/intel/mkl/8.0.1/lib/32/libmkl_lapack.a
BLAS=/opt/intel/mkl/8.0.1/lib/32/libmkl_ia32.a
LIBS=$(LAPACK) $(BLAS)  $(GUIDE)  -lpthread -lsvml
SYS=nag
DEFS= $(DEFS_CDF) $(DEFS_MPI)
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








