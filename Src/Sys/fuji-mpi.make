#
FC=fuji90
FC_ASIS=$(FC)
####FFLAGS= -O -Kfap -f 2006 -f 2004
FFLAGS= -O -f 2006 -f 2004 
LDFLAGS=-Wl,-Map linkmap
FFLAGS_DEBUG= -g -O0 -f 2006 -f 2004 -Hsu #-Hsu -ARy2
LIBS= -lfrtscalapack -l1upblas -l1utools -l.fuji.aux -lfrtredist \
      -lfrtblacsF77init_MPI -lfrtblacs_MPI \
      -lfrtblacsF77init_MPI -lfrtblacs_MPI -llapack -lblas \
       -l1umpich -lfujiarg
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



