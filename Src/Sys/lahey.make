#
FC=lf95
FC_ASIS=$(FC)
#
FFLAGS= -O --warn --quiet --tpp --ntrace
FFLAGS_DEBUG= -g -O0  -chk --trace
LDFLAGS=--staticlink
LIBS=
COMP_LIBS=linalg.a
SYS=bsd
DEFS=
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



