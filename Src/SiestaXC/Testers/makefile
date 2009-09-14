#!!@LICENSE
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
ARCH_MAKE_DEFAULT=../../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
# Path to source files
VPATH=../../../Src/SiestaXC/Testers
#
INCFLAGS=-I../../MPI
#
default: test1 test2 test3
#
libSiestaXC.a: 
	(cd .. ; $(MAKE) \
                    "VPATH_ROOT=../../Src" \
	            "INCFLAGS=-I../MPI" \
                    "FFLAGS=$(FFLAGS)" module \
	       ; cp -fp siestaxc.mod Testers )
#
test1: libSiestaXC.a test1.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test1 test1.o -lSiestaXC -L..
#
test2: libSiestaXC.a test2.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test2 test2.o -lSiestaXC -L..
#
test3: libSiestaXC.a test3.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test3 test3.o -lSiestaXC -L..
#
clean:
	rm -f *.o *.*d *.a test1 test2 test3

