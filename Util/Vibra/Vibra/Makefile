# Makefile for fcbuild
#
# Pablo Ordejon, August 1998
#
.SUFFIXES: .f .o .a 
#
default: what fcbuild vibrator
#
# System dependencies are dealt with by an include file. See Sys directory
# for examples and instructions. The variable ARCH should be in the 
# shell environment. 
#
ARCH?=unknown
#
include Sys/${ARCH}.make
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
what:
	@echo
	@echo "Compilation architecture to be used:  $(ARCH)"
	@echo "If that is not what you want, give the correct"
	@echo "value to the variable ARCH in your shell environment."
	@echo $(COMMENTS)
	@echo
#
#
#
BUILD-SRCS= fcbuild.f chkdim.f recoor.f
BUILD-OBJS=$(BUILD-SRCS:.f=.o)

VIBRA-SRCS= vibrator.f chkdim.f recoor.f hermdp.f paste.f klines.f outbands.f \
            reclat.f parse.f
VIBRA-OBJS=$(VIBRA-SRCS:.f=.o)

#
# Routines using vibra.h
#
fcbuild.o vibrator.o:  vibra.h
#
# Routines using fdf calls.
#
fcbuild.o klines.o outbands.o recoor.o vibrator.o: fdf/fdfdefs.h
#
FDF= fdf/fdf.o fdf/io.o
$(FDF):
	(cd fdf ; $(MAKE) "FC=$(FC)" "FFLAGS=$(FFLAGS)" fdf.o io.o)
#
fcbuild:$(BUILD-OBJS) $(FDF) 
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD-OBJS) $(FDF) $(LIBS)

vibrator : $(VIBRA-OBJS) $(FDF) 
	$(FC) -o vibrator \
	       $(LDFLAGS) $(VIBRA-OBJS) $(FDF) $(LIBS)
#
# Basic recipe to compile
#
.f.o:
	$(FC) -c $(FFLAGS) $<
#
clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f fcbuild vibrator *.o 
	(cd fdf ; make clean)



