# 
# This file is part of the SIESTA package.
#
# Copyright (c) Fundacion General Universidad Autonoma de Madrid:
# E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
# and J.M.Soler, 1996-2006.
# 
# Use of this software constitutes agreement with the full conditions
# given in the SIESTA license, as signed by all legitimate users.
#
# Makefile for Vibra package
#
# Pablo Ordejon, August 1998
# Modified by Alberto Garcia, February 2006
#
.SUFFIXES:
.SUFFIXES: .f .o .a 
#
default: what copy_fdf fcbuild vibrator
#
# You should have compiled Siesta in the main Src directory before
# attempting to compile the Vibra package.
#
include ../../../Src/arch.make
#
# (Note: The following is needed only if INFLAGS is not included
#  in the compilation rules in arch.make, but it is kept here in
#  case a crippled version of arch.make is used.)
#
# The fdf module file also needs to be present in the current directory
# Unfortunately the naming convention for the fdf module is compiler specific
# and so you might need to edit this. Should work for Intel and IBM compilers.
#
copy_fdf:
	cp -fp ../../../Src/fdf.mod .
#
# Uncomment the following line for debugging support
#
FFLAGS=$(FFLAGS_DEBUG)
#
what:   
	@echo
	@echo "Compilation architecture to be used:  $(SIESTA_ARCH)"
	@echo
#
#
BUILD-SRCS= fcbuild.f chkdim.f recoor.f io.f
BUILD-OBJS=$(BUILD-SRCS:.f=.o)

VIBRA-SRCS= vibrator.f chkdim.f recoor.f hermdp.f paste.f klines.f outbands.f \
            reclat.f parse.f io.f
VIBRA-OBJS=$(VIBRA-SRCS:.f=.o)

#
# Routines using vibra.h
#
fcbuild.o vibrator.o:  vibra.h
#
FDF=../../../Src/libfdf.a         
INCFLAGS=-I../../../Src           # To pick up fdf
#
fcbuild: $(BUILD-OBJS)
	$(FC) -o fcbuild \
	       $(LDFLAGS) $(BUILD-OBJS) $(FDF)

vibrator : $(VIBRA-OBJS) 
	$(FC) -o vibrator \
	       $(LDFLAGS) $(VIBRA-OBJS) $(FDF)
#
clean: 
	@echo "==> Cleaning object and executable files"
	rm -f fcbuild vibrator *.o 



