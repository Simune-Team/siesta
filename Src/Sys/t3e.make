#
# Makefile include file for Cray's T3E
# Alberto Garcia <wdpgaara@lg.ehu.es>, Aug 2, 1999
#
# Compiler invocation. Note '-em' to produce module files and
# '-dp' to disable double precision
#
FC=f90 -em -dp
#
# To compile some interface modules, we need the actual meaning
# of 'double precision' on the Cray.
#
FC_ASIS=f90 -em
#
# Whatever needed, except -dp and -em... 
#
FFLAGS=  -O scalar2,pipeline2,aggress
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -Rabc -ei
#
#
#
LIBS= -lsci -lmpi 
#LIBS= -lsci -lmpi -lapp
#
SYS=t3e
RANLIB=
#
# Location of mpif.h include file
#
MPI_INCLUDE=/usr/local/include
#
# Definition to trigger MPI conditional compilation
#
DEFS=-DMPI -DCRAY
MPILIB=libmpi_f90.a
#
# Actual compilation recipes for siesta code.
# Specify "-p ." to let the compiler know that
# the modules are in the current directory.
# It could work without it...
#
.F.o:
	$(FC) -c $(FFLAGS) -p . $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) -p .  $<
.F90.o:
	$(FC) -c $(FFLAGS) -p .  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) -p .  $<
#
