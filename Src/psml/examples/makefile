.SUFFIXES: 
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default:  test_psml
#
ARCH_MAKE_DEFAULT=../src/fortran.mk
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)
#
OBJFILES= test_psml.o psml_die.o m_getopts.o\
          interpolate_drh.o
#
# Uncomment the following if you implement a new interpolator
#
#EXTRA_OBJS=interpolate_other.o
#
INCFLAGS= -I../src
LIBS= ../src/libpsml.a
#
test_psml: $(OBJFILES) $(EXTRA_OBJS) $(LIBS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o test_psml $(INCFLAGS) \
              $(OBJFILES) $(EXTRA_OBJS) $(LIBS)
#
test_psml.o: m_getopts.o
#
clean:
	rm -f *.o test_psml
