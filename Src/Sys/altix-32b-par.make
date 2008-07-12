SIESTA_ARCH=IFORT 32bits PARALLEL
#
FC= ifort

FFLAGS_DEBUG= -g
FFLAGS=-O3 -cpp              # -assume nounderscore
FFLAGS=-g -cpp -warn all -WB -C -debug full -traceback -fpe0
LDFLAGS= 
COMP_LIBS=
RANLIB=echo

#
MPI_INTERFACE=    libmpi_f90.a
MPI_INCLUDE=.
MPI_LIBS=         #-lblacsgm
DEFS_MPI=         -DMPI



#BLAS = -L/gpfs/apps/LAPACK/lib32 -llapack \
#	     -L/gpfs/apps/SCALAPACK/lib32 -lblas


#BLAS= -L/gpfs/apps/SCALAPACK/lib32 -lscalapack \
#        /gpfs/apps/SCALAPACK/lib32/blacsF77init_MPI-PPC-0.a \
#        /gpfs/apps/SCALAPACK/lib32/blacsCinit_MPI-PPC-0.a \
#        /gpfs/apps/SCALAPACK/lib32/blacs_MPI-PPC-0.a \
#      -L/gpfs/apps/LAPACK/lib32 -llapack \
#      -L/gpfs/apps/SCALAPACK/lib32 -lblas

BLAS=-lscs
BLAS= /apps/LAPACK/3.1.1/lib/liblapack.a                   \
      /apps/SCALAPACK/1.8.0/lib/libscalapack.a             \
      /apps/SCALAPACK/1.8.0/lib/blacsF77init_MPI-LINUX-0.a \
      /apps/SCALAPACK/1.8.0/lib/blacsCinit_MPI-LINUX-0.a   \
      /apps/SCALAPACK/1.8.0/lib/blacs_MPI-LINUX-0.a        \
      /apps/SCALAPACK/1.8.0/lib/blas_LINUX.a
#PAPI=-L/gpfs/apps/PAPI/papi-3.2.1/32/lib/ -lpapi
LIBS= $(PAPI) $(BLAS) $(MPITRACER) -lmpi -lgfortran

SYS=xlf
DEFS= $(DEFS_MPI) $(DEFS_CDF) $(DEFS_TRACE) -DDEBUG -DALTIX
FREE_F90=-qsuffix=f=f90
#
.F90.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS)   $<
.F.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
#.f.o:
#	$(FC) -qsuffix=f=f -c $(INCFLAGS) $(FFLAGS)   $<
#
.f.o:
	$(FC) -c $(INCFLAGS) $(FFLAGS) $(DEFS) $<
