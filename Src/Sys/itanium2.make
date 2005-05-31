SIESTA_ARCH=itanium2
#
# Serial compilation
#
FC=efc -Vaxlib
#
FC_ASIS=$(FC)
#
FFLAGS= -O2 -tpp2 -W0
FFLAGS_DEBUG= -g -O0
LDFLAG=#-Vaxlib 
RANLIB=echo
LAPACK=-L/opt/intel/mkl72/lib/64 -lmkl_lapack64 -lmkl -lguide -lpthread
LIBS=$(LAPACK)
SYS=bsd
DEFS= -DWXML_INIT_FIX -DALLOC_SAVE_BUG
COMP_LIBS=
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



