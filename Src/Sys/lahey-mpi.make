#
FC=lf95
FC_ASIS=$(FC)
FFLAGS= -O  --ntrace --tpp
LDFLAGS=
FFLAGS_DEBUG= -g -O0 -chk
LIBS= -lscalapack-lf95 -lpblas-lf95 -ltools-lf95 -lredist-lf95 \
      -lblacsF77init_MPI-lf95-0  -lblacs_MPI-lf95-0 \
      -lblacsF77init_MPI-lf95-0  -lblacs_MPI-lf95-0 \
      -llapack-lf95 -lblas-lf95 \
      -lmpich-lf95
SYS=bsd
MPILIB=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
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



