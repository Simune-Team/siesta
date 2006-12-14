#
# Prototype .make file for Mare Nostrum in Barcelona
# (Courtesy of Manuel Cobian <cobian@icmab.es>)
# Some tuning might still be needed.
#
SIESTA_ARCH=XLF 32bits PARALLEL
#
FC=mpif90 #xlf90_r
#
FFLAGS_DEBUG= -g
FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q32
FFLAGS_parse=-qsuffix=f=f -qfixed
LDFLAGS= -q32
COMP_LIBS=           # Might need fortran D&C routines
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=
MPI_LIBS= #-lblacsgm
DEFS_MPI= -WF,-DMPI
#

#
BLAS= -L/gpfs/apps/SCALAPACK/lib32 -lscalapack \
        /gpfs/apps/SCALAPACK/BLACS/BLACS/LIB/blacsF77init_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/BLACS/BLACS/LIB/blacsCinit_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/BLACS/BLACS/LIB/blacs_MPI-PPC-0.a \
      -L/gpfs/apps/LAPACK/lib32 -llapack \
      -L/gpfs/apps/SCALAPACK/BLAS/lib32 -lblas
##MPITRACER= -L/gpfs/apps/CEPBATOOLS/lib/32 -lmpitrace
#
##LIBS=  $(MPITRACER) $(BLAS)
LIBS=  $(BLAS)
SYS=xlf
DEFS= $(DEFS_MPI) $(DEFS_CDF)
FREE_F90=-qsuffix=f=f90
#
.F90.o:
	$(FC) -qsuffix=cpp=F90 -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -qsuffix=f=f90 -c $(FFLAGS) $(INCFLAGS)   $<
.F.o:
	$(FC) -qsuffix=cpp=F -qfixed -c $(FFLAGS) $(INCFLAGS) $(DEFS) $<
.f.o:
	$(FC) -qsuffix=f=f -qfixed -c $(FFLAGS) $(INCFLAGS)   $<
#
