#
FC=pgf90
FC_ASIS=$(FC)
#
FFLAGS=  -fast
FFLAGS_DEBUG= -g -O0
LIBS= -L/usr/local/lib -llapack -lblas -lg2c
SYS=bsd
#
MPILIB=
MPI_INCLUDE=/usr/local/include
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



