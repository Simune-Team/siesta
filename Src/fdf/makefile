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
libfdf.a:	fdf.o fdf_mod.o
	ar qv libfdf.a fdf.o fdf_mod.o
fdf_mod.o: fdf_mod.f
	$(FC) $(FFLAGS) -c fdf_mod.f
	@cp -p *.mod ..
fdf.o:  fdf.f 
	$(FC) $(FFLAGS) -c fdf.f
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



