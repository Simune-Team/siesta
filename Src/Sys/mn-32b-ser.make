SIESTA_ARCH=XLF 32bits SERIAL
#
FC=xlf90   #xlf90_r
#
FFLAGS_DEBUG= -g
#FFLAGS=-O3 -qstrict -qtune=ppc970 -qarch=ppc970 -q32
FFLAGS=-g -qfullpath -qtune=ppc970 -qarch=ppc970 -q32
FFLAGS_parse=-qsuffix=f=f -qfree #-qfixed
LDFLAGS= -q32
#LDFLAGS= -q32 -g
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
#MPI_INTERFACE=    libmpi_f90.a
#MPI_INCLUDE=.
#MPI_LIBS=         #-lblacsgm
#DEFS_MPI=         -WF,-DMPI



#BLAS = -L/gpfs/apps/LAPACK/lib32 -llapack \
#	     -L/gpfs/apps/SCALAPACK/lib32 -lblas


BLAS= -L/gpfs/apps/SCALAPACK/lib32 -lscalapack \
        /gpfs/apps/SCALAPACK/lib32/blacsF77init_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib32/blacsCinit_MPI-PPC-0.a \
        /gpfs/apps/SCALAPACK/lib32/blacs_MPI-PPC-0.a \
      -L/gpfs/apps/LAPACK/lib32 -llapack \
      -L/gpfs/apps/SCALAPACK/lib32 -lblas

#TRACE_DIR=/gpfs/apps/CEPBATOOLS/mpitrace/32/lib/
#MPITRACER=-L$(TRACE_DIR) -lmpitracef -lxml2
#  -lmpitrace
#DEFS_TRACE= -WF,-DMPI_TRACE

PAPI=-L/gpfs/apps/PAPI/papi-3.2.1/32/lib/ -lpapi
LIBS= $(PAPI) $(BLAS) $(MPITRACER)
#LIBS= $(BLAS)

SYS=xlf
DEFS_DEBUG= -WF,-DDEBUG
DEFS= $(DEFS_MPI) $(DEFS_CDF) $(DEFS_TRACE) $(DEFS_DEBUG) -WF,-DXLF
FREE_F90=-qsuffix=f=f90
#
.F90.o:
	$(FC) -qsuffix=cpp=F90 -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -qsuffix=f=f90 -c $(INCFLAGS) $(FFLAGS)   $<
.F.o:
	$(FC) -qsuffix=cpp=F -c $(INCFLAGS) -qfixed $(FFLAGS) $(DEFS) $<
#.f.o:
#	$(FC) -qsuffix=f=f -qfixed -c $(INCFLAGS) $(FFLAGS)   $<
#
.f.o:
	$(FC) -qsuffix=cpp=f -c $(INCFLAGS) -qfixed $(FFLAGS) $(DEFS) $<
