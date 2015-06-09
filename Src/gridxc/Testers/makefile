#!!@LICENSE
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90
#
ARCH_MAKE_DEFAULT=../lib-arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
VPATH=.
#
GRIDXC_ROOT=..
GRIDXC_LIBS=$(GRIDXC_ROOT)
GRIDXC_INCLUDE=$(GRIDXC_ROOT)
#
INCFLAGS=-I$(GRIDXC_INCLUDE)  
#
default: test1 test2 test3 test4
#
LIBS = -L $(GRIDXC_LIBS) -lGridXC
#
test1:  test1.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test1 test1.o handlers.o $(LIBS)
#
test2:  test2.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test2 test2.o handlers.o $(LIBS)
#
test3:  test3.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test3 test3.o handlers.o $(LIBS)
#
test4:  test4.o handlers.o
	$(FC) $(FFLAGS) $(DEFS) $(INCFLAGS) -o test4 test4.o handlers.o $(LIBS)
#
clean:
	rm -f *.o *.*d *.a test1 test2 test3 test4
	rm -f atomVxc.out cellVxc.out  debugXC*out
	rm -f vdw_kernel.table*  *.Vxc


