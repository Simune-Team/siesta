#
#  Makefile for fdf example.
#
FC?=f77
CFLAGS= 
FFLAGS=       $(CFLAGS)
LDFLAGS=      $(CFLAGS)
#
OBJS=	sample.o fdf.o io.o
#
sample: $(OBJS)
	$(FC) $(LDFLAGS) -o sample $(OBJS)
#
check:	sample sample.fdf
	sample < sample.fdf
#
io_sample: io_sample.o io.o
	$(FC) $(LDFLAGS) -o io_sample io_sample.o io.o
#
io_check: io_sample
	io_sample
#
module:  libfdf.a
libfdf.a:	fdf.o 
	ar qv libfdf.a fdf.o
fdf.o:  fdf.f fdf_mod.f
	$(FC) $(FFLAGS) -c fdf_mod.f
	$(FC) $(FFLAGS) -c fdf.f
	@cp -p *.mod ..
#
clean:
	rm -f $(OBJS) libfdf.a *.mod fdf_mod.o io_sample.o sample io_sample
	rm -f OUTPUT FDF.debug file2 file3 sample.out
#
.f.o:
	$(FC) $(FFLAGS) -c $<
#
fdf.o:		fdf.h
sample.o:	fdfdefs.h
#



