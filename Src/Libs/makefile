#
#  Makefile for Libraries
#  Some systems only need the stuff in dc_lapack.f (Divide and Conquer)
#
default: linalg.a
#
include ../arch.make
#
ALL_OBJS= lapack.o dc_lapack.o blas.o machine.o
# Add machine.o here only if needed...
DC_OBJS=  dc_lapack.o   # machine.o
#
default: linalg.a
#
libblas.a: blas.o
	ar $(ARFLAGS_EXTRA) cru libblas.a blas.o
	$(RANLIB) libblas.a
	mv libblas.a ..
liblapack.a: lapack.o machine.o
	ar $(ARFLAGS_EXTRA) cru liblapack.a lapack.o machine.o
	$(RANLIB) liblapack.a
	mv liblapack.a ..
linalg.a: $(ALL_OBJS)
	ar $(ARFLAGS_EXTRA) cru linalg.a $(ALL_OBJS)
	$(RANLIB) linalg.a
	mv linalg.a ..
dc_lapack.a: $(DC_OBJS)
	ar $(ARFLAGS_EXTRA) cru dc_lapack.a $(DC_OBJS)
	$(RANLIB) dc_lapack.a
	mv dc_lapack.a ..
#
clean:
	rm -f $(ALL_OBJS) *.a




