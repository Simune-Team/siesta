.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default:  test_ncps
#
ARCH_MAKE_DEFAULT=../src/fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
VPATH=.
#
OBJFILES= test_ncps.o m_getopts.o die.o psml_die.o
INCFLAGS= -I../src
LIBS= ../src/libncps.a $(PSML_LIBS) 
#
test_ncps: $(OBJFILES) $(LIBS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o test_ncps $(INCFLAGS) \
              $(OBJFILES) $(LIBS)
#
test_ncps.o: m_getopts.o
#
clean:
	rm -f *.o test_ncps
