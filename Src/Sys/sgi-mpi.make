#
FC=f90 -n32
#
FFLAGS=  -O3 
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -O0
#
LIBS= -lcomplib.sgimath libmpiblacsn32.a libscalapackn32.a libmpiblacsn32.a -lmpi 
SYS=bsd
#
# Location of mpif.h include file
#
MPILIB=libmpi_f90.a
MPI_INCLUDE=/usr/include
#
# Definition to trigger MPI conditional compilation
#
DEFS=-DMPI
#
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<
#
