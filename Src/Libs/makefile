#
#  Makefile for Libraries
#  Some systems only need the stuff in dc_lapack.f (Divide and Conquer)
#
include ../arch.make
#
FC?=f77
CFLAGS?= 
FFLAGS?=       $(CFLAGS)
LDFLAGS?=      $(CFLAGS)
#
ALL_OBJS= lapack.o dc_lapack.o blas.o machine.o
DC_OBJS=  dc_lapack.o
#
default: linalg.a
#
# DO NOT OPTIMIZE THIS!
#
machine.o:
	$(FC) -c machine.f
#
linalg.a:  $(ALL_OBJS)
	ar qv linalg.a $(ALL_OBJS)
	$(RANLIB) linalg.a
	mv linalg.a ..
dc_lapack.a:  $(DC_OBJS)
	ar qv dc_lapack.a $(DC_OBJS)
	$(RANLIB) dc_lapack.a
	mv dc_lapack.a ..
#
clean:
	rm -f $(ALL_OBJS) *.a




