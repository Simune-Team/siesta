#
FC=ff90
FC_ASIS=$(FC)
#
FFLAGS= -fast
FFLAGS_DEBUG= -g -O0
LIBS= -lfrtscalapack -l1upblas -l1utools -l.pgi.aux -lfrtredist \
      -lfrtblacsF77init_MPI -lfrtblacs_MPI \
      -lfrtblacsF77init_MPI -lfrtblacs_MPI -llapack -lblas \
       -l1umpich -lpgiarg
SYS=bsd
#
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





