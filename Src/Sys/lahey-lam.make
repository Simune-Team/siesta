#
FC=lf95
FC_ASIS=$(FC)
FFLAGS= -O  --ntrace --tpp
LDFLAGS=--staticlink
FFLAGS_DEBUG= -g -O0
LIBS= -lscalapack_lf95-lam -lpblas_lf95-lam -ltools_lf95-lam -lredist_lf95-lam\
      -lblacsF77init_MPI-lf95-lam-0  -lblacs_MPI-lf95-lam-0 \
      -lblacsF77init_MPI-lf95-lam-0  -lblacs_MPI-lf95-lam-0 \
      -llapack-lf95 -lblas-lf95 \
      -L/usr/local/lam/lib -lmpi -ltstdio -ltrillium -largs -lt
SYS=bsd
MPILIB=libmpi_f90.a
MPI_INCLUDE=/usr/local/lam/include
DEFS=-DMPI
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



