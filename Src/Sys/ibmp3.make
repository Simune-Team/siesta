SIESTA_ARCH=ibmp3
#
# IBM without ESSL
#
FC=xlf -bloadmap:MAP
FC_ASIS=$(FC)
#
FFLAGS= -O3 -qmaxmem=16000 -qtkq_size=25000 -qnolm
FFLAGS_DEBUG= -g -qmaxmem=16000 -qtkq_size=25000 -qnolm
# -qipa gives speed optimization of a few percent, but takes long to link
LIBS= -Wl,-bD:2000000000 -qnolm
COMP_LIBS=linalg.a
MPILIB=
SYS=ibm
RANLIB=ranlib
MPI_INCLUDE=/usr/local/include
DEFS=
#
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<
