#
#  Makefile for Libraries
#
FC?=f77
CFLAGS= 
FFLAGS=       $(CFLAGS)
LDFLAGS=      $(CFLAGS)
#
default: linalg.a
OBJS=	lapack.o blas.o
#
linalg.a:  $(OBJS)
	ar qv linalg.a $(OBJS)
	@mv linalg.a ..
#
clean:
	rm -f $(OBJS) *.a
#
.f.o:
	$(FC) $(FFLAGS) -c $<
#




